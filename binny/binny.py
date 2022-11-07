#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 22 10:50:35 2021

@author: oskar.hickl
"""
import logging
import sys
import os
import glob
from binny.utils import *


def binny():


    sm_pur = snakemake.params['purity']
    sm_m_comp = snakemake.params['min_completeness']
    sm_s_comp = snakemake.params['start_completeness']
    sm_mask = snakemake.params['mask_disruptive_sequences']
    sm_scmag = snakemake.params['extract_scmags']
    sm_nx = snakemake.params['nx_val']
    sm_min_co = snakemake.params['min_cutoff']
    sm_max_co = snakemake.params['max_cutoff']
    sm_min_m_co = snakemake.params['min_cutoff_marker']
    sm_max_m_co = snakemake.params['max_cutoff_marker']
    sm_n_cont = snakemake.params['max_n_contigs']
    sm_mark_lvl = snakemake.params['max_marker_lineage_depth_lvl']
    sm_n_tries = snakemake.params['max_embedding_tries']
    sm_i_depth = snakemake.params['include_depth_initial']
    sm_m_depth = snakemake.params['include_depth_main']
    sm_e_range = snakemake.params['hdbscan_epsilon_range']
    sm_m_s_range = snakemake.params['hdbscan_min_samples_range']
    sm_write = snakemake.params['write_contig_data']

    binny_out = snakemake.params['binny_out']
    sample = snakemake.params['sample']
    mg_depth_file = snakemake.input['mgdepth']
    assembly = snakemake.input['assembly']
    annot_file = snakemake.params['gff']
    raw_annot = snakemake.input['raw_gff']
    tigrfam2pfam_file = snakemake.params['t2p']
    taxon_marker_set_file = snakemake.params['marker_sets']
    prokka_checkm_marker_hmm_out = snakemake.input['hmm_markers']
    functions = snakemake.params['py_functions']
    min_purity = float(sm_pur) if type(sm_pur) == str else sm_pur
    min_completeness = float(sm_m_comp) if type(sm_m_comp) == str else sm_m_comp
    starting_completeness = float(sm_s_comp) if type(sm_s_comp) == str else sm_s_comp
    kmers = snakemake.params['kmers']
    mask_disruptive_sequences = eval(sm_mask) if type(sm_mask) == str else sm_mask
    extract_scmags = eval(sm_scmag) if type(sm_scmag) == str else sm_scmag
    coassembly_mode = snakemake.params['coassembly_mode']
    nx_val = int(sm_nx) if type(sm_nx) == str else sm_nx
    min_contig_length = int(sm_min_co) if type(sm_min_co) == str else sm_min_co
    max_contig_length = int(sm_max_co) if type(sm_max_co) == str else sm_max_co
    min_contig_length_marker = int(sm_min_m_co) if type(sm_min_m_co) == str else sm_min_m_co
    max_contig_length_marker = int(sm_max_m_co) if type(sm_max_m_co) == str else sm_max_m_co
    max_contig_threshold = float(sm_n_cont) if type(sm_n_cont) == str else sm_n_cont
    max_marker_lineage_depth_lvl = int(sm_mark_lvl) if type(sm_mark_lvl) == str else sm_mark_lvl
    max_embedding_tries = int(sm_n_tries) if type(sm_n_tries) == str else sm_n_tries
    include_depth_initial = eval(sm_i_depth) if type(sm_i_depth) == str else sm_i_depth
    include_depth_main = eval(sm_m_depth) if type(sm_m_depth) == str else sm_m_depth
    hdbscan_epsilon_range = [float(epsilon) for epsilon in sm_e_range.split(',')] if type(sm_e_range) == str else sm_e_range
    hdbscan_min_samples_range = [int(min_sample) for min_sample in sm_m_s_range.split(',')] if type(sm_m_s_range) == str else sm_m_s_range
    dist_metric = snakemake.params['distance_metric']
    write_contig_data = eval(sm_write) if type(sm_write) == str else sm_write

    intermediary_file_dir = 'intermediary'

    threads = snakemake.threads
    log = snakemake.log[0]

    n_dim = 2

    sys.path.append(functions)

    # To achieve reproducible results with HDBSCAN and ensure same seed, because other tools that accept seed arguments,
    # might mess with the global numpy seed
    np.random.seed(0)

    logging.basicConfig(filename=log, level=logging.INFO, format='%(asctime)s - %(message)s',  # logging.INFO
                        datefmt='%d/%m/%Y %I:%M:%S %p', filemode='w')

    numba_logger = logging.getLogger('numba')
    numba_logger.setLevel(logging.WARNING)

    logging.info('Starting Binny run for sample {0}.'.format(sample))


    # Check if bin dir empty
    if glob.glob(os.path.join(binny_out, "bins/*.fasta")):
        logging.error('Bin dir contains fasta files. Move or delete. Exiting.')
        raise Exception

    # Load TIGRFAMs to PFAMs conversion table.
    tigrfam2pfam_data = tigrfam2pfam_dict(tigrfam2pfam_file)

    # Merge Prokka gff with marker set data, load annotation df, and load assembly.
    checkm_hmmer_search2prokka_gff(prokka_checkm_marker_hmm_out, raw_annot, gff_out_path=os.path.join(binny_out,
                                                                                                    'intermediary'))
    annot_df, annot_dict = gff2ess_gene_df(annot_file, target_attribute='checkm_marker', get_dict=True)
    assembly_dict = load_fasta(assembly)

    # Build marker set graph db.
    taxon_marker_sets = load_checkm_markers(taxon_marker_set_file)

    # Look for complete genomes on single contigs
    all_good_bins = {}
    if extract_scmags:
        logging.info('Looking for single contig bins.')
        single_contig_bins = get_single_contig_bins(annot_df, all_good_bins, n_dim, taxon_marker_sets, tigrfam2pfam_data,
                                                    threads, max_marker_lineage_depth_lvl=max_marker_lineage_depth_lvl)
        logging.info('Found {0} single contig bins.'.format(len(single_contig_bins)))
    else:
        single_contig_bins = []


    logging.info(f'Calculating N{nx_val}'.format(len(single_contig_bins)))
    nx = calc_assembly_nx(assembly_dict, single_contig_bins, nx_val)
    nx2 = calc_assembly_nx(assembly_dict, [], nx_val)
    logging.info(f'N{nx_val} is {nx}, with scMAGs would be {nx2}.'.format(len(single_contig_bins)))
    min_contig_length = min(max(nx, min_contig_length), max_contig_length)
    min_contig_length_marker = min(max(int(nx / 3), min_contig_length_marker), max_contig_length_marker)


    contig_list = [[contig] + [seq] for contig, seq in assembly_dict.items() if contig not in single_contig_bins
                                                                                and len(seq) >= 500]

    contig_rrna_crispr_region_dict = gff2low_comp_feature_dict(annot_file)
    if mask_disruptive_sequences:
        logging.info('Masking potentially disruptive sequences from k-mer counting.')
        mask_rep_featrues(contig_rrna_crispr_region_dict, contig_list)
    else:
        logging.info('Not masking potentially disruptive sequences from k-mer counting.')

    # Get length normalized k-mer frequencies.
    kmer_sizes = [int(kmer) for kmer in kmers.split(',')]

    start = timer()
    kfreq_array = get_contig_kmer_matrix(contig_list, kmer_sizes, threads)
    end = timer()
    logging.info('K-mer frequency matrix created in {0}s.'.format(int(end - start)))

    # Make array, removing fully masked sequences with no counts and standardize k-mer freq data
    x = np.array([c_kfreq[1:] for c_kfreq in kfreq_array[1:] if not sum(c_kfreq[1:]) == 0])
    x_contigs = [c_kfreq[0] for c_kfreq in kfreq_array[1:] if not sum(c_kfreq[1:]) == 0]

    main_contig_data_dict = {cont: seq for cont, seq in zip(x_contigs, x)}

    # Load depth data
    depth_dict = load_depth_dict(mg_depth_file)
    n_depth_samples = list(depth_dict.values())[0].shape[0]
    logging.debug(f'list(depth_dict.values())[0].shape[0]: {n_depth_samples}')

    if coassembly_mode == 'on' or (coassembly_mode == 'auto' and n_depth_samples > 1):
        min_contig_length_marker = 500
        if n_depth_samples >= 3:
            include_depth_main = True
        logging.info('Using coassembly mode.')
    else:
        logging.info('Using single sample mode.')

    # Run iterative dimension reduction, manifold learning, cluster detection and assessment.
    all_good_bins, contig_data_df_org, min_purity = iterative_embedding(x_contigs, depth_dict, all_good_bins, starting_completeness,
                                                            min_purity, min_completeness, threads, n_dim, annot_file,
                                                            mg_depth_file, single_contig_bins, taxon_marker_sets,
                                                            tigrfam2pfam_data, main_contig_data_dict, assembly_dict,
                                                            max_contig_threshold, min_contig_length_marker,
                                                            include_depth_initial, max_embedding_tries,
                                                            include_depth_main, hdbscan_epsilon_range,
                                                            hdbscan_min_samples_range, dist_metric,
                                                            contigs2clusters_out_path=os.path.join(binny_out, 'intermediary'),
                                                            max_marker_lineage_depth_lvl=max_marker_lineage_depth_lvl,
                                                            coassembly_mode=coassembly_mode)

    all_contigs = []
    for bin in all_good_bins:
        all_contigs.extend(all_good_bins[bin]['contigs'])
    if len(all_contigs) != len(set(all_contigs)):
        logging.warning('WARNING: {0} duplicate contigs in bins found!'.format(len(all_contigs) - len(set(all_contigs))))

    # Write bin fastas.
    write_bins(all_good_bins, assembly, min_comp=int(min_completeness), min_pur=int(min_purity),
            bin_dir=os.path.join(binny_out, 'bins'))

    bin_dict = {contig: bin for bin, bin_data in all_good_bins.items() for contig in bin_data['contigs']}

    all_cont_data_dict = {}

    for contig, k_freq in main_contig_data_dict.items():
        all_cont_data_dict[contig] = {'bin': bin_dict.get(contig, 'N'),
                                    'k-mer_freqs': ';'.join([str(k) for k in list(k_freq)]),
                                    'depths': ';'.join([str(d) for d in list(depth_dict.get(contig))])}

    all_cont_data_dict.update({contig: {'bin': contig,
                                    'k-mer_freqs': '',
                                    'depths': ';'.join([str(d) for d in list(depth_dict.get(contig))])}
                            for contig in single_contig_bins})

    contig_data_df = pd.DataFrame.from_dict(all_cont_data_dict, orient='index', columns=['bin', 'k-mer_freqs', 'depths'])
    # Downcast df
    downcast_pd_df(contig_data_df)

    if write_contig_data:
        logging.info('Writing contig data to file.')
        compression_opts = dict(method='gzip')
        contig_data_df.to_csv(os.path.join(binny_out, 'contig_data.tsv.gz'), header=True, index=True, index_label='contig',
                            chunksize=250000, compression=compression_opts, sep='\t')

    os.mknod(os.path.join(binny_out, "binny.done"))

    logging.info('Run finished.')
