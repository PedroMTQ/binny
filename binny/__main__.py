#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Based on https://github.com/my-yangj/2019-snakemake-cli/blob/master/cli/command.py
Command line interface driver for snakemake workflows
"""
import argparse
import os.path
import sys
import tarfile
import pandas as pd
import snakemake
from distutils.sysconfig import get_python_lib

from binny.utils.utils import run_cmd

try:
    from mantis.assembler import check_installation, setup_databases
    from mantis.utils import download_file
    from mantis.exceptions import InstallationCheckNotPassed
except (ImportError, ModuleNotFoundError):
    env_path = '/'.join(sys.executable.split('/')[:-2])
    python_version = '.'.join([str(sys.version_info.major), str(sys.version_info.minor)])
    sys.path.append(os.path.join(f'{env_path}/lib/python{python_version}/site-packages/mantis'))
    from assembler import check_installation, setup_databases
    from exceptions import InstallationCheckNotPassed
    from utils import download_file


BINNY_DIR = os.path.abspath(os.path.dirname(__file__))
snakefile = os.path.join(BINNY_DIR, "Snakefile")
def_configfile = os.path.join(BINNY_DIR, 'config/config.default.yaml')
int_config = os.path.join(BINNY_DIR, 'config/binny_snakemake_defaults.yaml')
src_dir = os.path.join(BINNY_DIR, 'workflow/scripts')
bin_dir = os.path.join(BINNY_DIR, 'workflow/bin')
config_dir = os.path.join(BINNY_DIR, 'config')
conda_dir = os.path.join(BINNY_DIR, 'conda')
cwd = os.getcwd()

if src_dir not in sys.path:
    sys.path.append(src_dir)
    import binny.utils.remove_unused_checkm_hmm_profiles as prepCheckM


def arg_parser():
    parser = argparse.ArgumentParser(description='Run binny on target assembly.',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # Either use config file only or specify manually. Required.
    input_type = parser.add_mutually_exclusive_group(required=True)

    input_type.add_argument('-a', '--assembly', type=str, help='Path to an assembly fasta.')

    # Required coverage data.
    cov_args = parser.add_mutually_exclusive_group(required=True)
    cov_args.add_argument('-b', '--bam', type=str, nargs='*',
                          help='Path to a bam file(s) to calculate depth from. Use wildcards for multiple samples, '
                               'e.g.: "path/to/my/mappings/*.bam" '
                               'or "path/to/my/mappings/with/*/different/folder/*/structure/*.bam" '
                               'or "path/to/my/mappings/my_mapping.bam. Leave empty if you have an average depth per '
                               'contig file to supply to binny.')
    cov_args.add_argument('--contig_depth', type=str, default=None,
                          help='Path to an average depth per contig tsv file. Leave empty if you supply a bam file for '
                               'binny to calculate average contig depth from. First column needs to be the contig ids, '
                               'subsequent column(s) for depht(s).')

    # Remaining required input.
    out_args = parser.add_argument_group('Output', 'Other arguments required to run binny.')
    out_args.add_argument('-o', '--outputdir', type=str, required=True,
                          help='Path to desired output dir binny should create and store results in.')
    out_args.add_argument('--sample', type=str, default='sample', help='Sample name.')
    out_args.add_argument('--tmp_dir', type=str, default=None,
                          help='Path to a temporary directory to write to. Defaults to outputdir/tmp')

    # Resource params.
    res_args = parser.add_argument_group('Computing resources', 'Parameters concerning computing resoures.')
    res_args.add_argument('-t', '--threads', type=int, default=1,
                          help='Maximum number of cpus/threads to use.')
    res_args.add_argument('--big_mem_avail', action='store_true', default=False,
                          help='Set to use high memory capacity computer node.')
    res_args.add_argument('--big_mem_per_core_gb', type=int, default=26,
                          help='Specify the amount of memory per core for high memory capacity node. '
                               'Use with \'--big-mem\'.')
    res_args.add_argument('--normal_mem_per_core_gb', type=int, default=2,
                          help='Specify the amount of memory per core for your system.')



    # Snakemake arguments
    sm_args = parser.add_argument_group('Snakemake arguments', 'Arguments concerning the execution of Snakemake.')
    sm_args.add_argument('-c', '--config_file', type=str, default=None,
                         help=f'Path to config file. A default config file to use as template '
                              f'can be found at: {def_configfile}.')
    sm_args.add_argument('-m', '--use_cluster', action='store_true', default=False,
                         help='Use cluster to submit jobs to instead of running locally.')
    sm_args.add_argument('-csc', '--cluster_submission_command', type=str,
                         default='{cluster.call} {cluster.runtime}{resources.runtime} '
                                 '{cluster.mem_per_cpu}{resources.mem} {cluster.nodes} {cluster.qos} '
                                 '{cluster.threads}{threads} {cluster.partition} {cluster.stdout}',
                         help='Submission command of a cluster or batch system to use (Atm only slurm support.')
    sm_args.add_argument('-spp', '--scheduler_preset_path', type=str, default=os.path.join(config_dir,
                                                                                           'slurm.config.yaml'),
                         help='Path to the scheduler preset file for cluster submission.')
    sm_args.add_argument('-jn', '--job_name', type=str, default='binny.{rulename}.{jobid}.sh',
                         help='Naming scheme for cluster job scripts.')
    sm_args.add_argument('-r', '--report', action='store_true', default=False,
                         help='Generate a report (it\'s recommended to run -c, -f and -l with -r).')
    sm_args.add_argument('-d', '--dry-run', action='store_true', default=False, help='Perform a dry-run.')
    sm_args.add_argument('-u', '--unlock', action='store_true', default=False,
                         help='Unlock working directory (only necessary for crash/kill recovery).')
    sm_args.add_argument('-f', '--force', action='store_true', default=False,
                         help='Force all output files to be re-created.')

    # Other resource arguments
    env_args = parser.add_argument_group('Environment arguments', 'Arguments concerning the execution environment.')
    env_args.add_argument('--db_path', type=str, default=os.path.join(BINNY_DIR, 'database'),
                          help='Absolute path to put binny dbs in. If left empty they will be put into \'database\' in '
                               'the binny main dir.')
    env_args.add_argument('--conda_env_path', type=str, default=os.path.join(BINNY_DIR, 'conda'),
                          help='Set path for conda env to be installed to. By default, they will be put in `conda` '
                               'in the binny dir.')

    # binny parameters
    bi_args = parser.add_argument_group('binny arguments', 'Arguments concerning parameters for binny.')
    bi_args.add_argument('--kmers', type=str, default='2,3,4', help='Input a list of kmers, e.g. \'2,3,4\'.')
    bi_args.add_argument('--mask_disruptive_sequences', action='store_true', default=True,
                         help='Toggle masking of potentially disruptive contig regions (e.g. rRNA and CRISPR elements) '
                              'from k-mer counting')
    bi_args.add_argument('--extract_scmags', action='store_true', default=True,
                         help='Extract single contig MAGs of at least 90%% purity and 92.5%% completeness.')

    bi_args.add_argument('--coassembly_mode', type=str, default='auto', choices=['auto', 'on' , 'off'],
                         help='Will use coassembly mode, starting with contigs >= 500 bp instead of high threshold, '
                              'decreasing, if set to \'on\' or if \'auto\' and multiple depth files are detected. '
                              'Choose between, \'auto\', \'on\' , \'off\'.')
    bi_args.add_argument('--NX_value', type=int, default=90,
                         help='binny prefilters assemblies based on N<X> value to try and take as much information as '
                              'possible into account,  while minimizing the amount of noise. Be aware that, depending '
                              'on the assembly quality, low values as the N<X> might results in leaving out a large '
                              'proportion of the assembly (if the max_cont_length cutoffs are set high as well).')
    bi_args.add_argument('--min_cont_length_cutoff', type=int, default=2250,
                         help='Minimum contig length. Caps NX filtering value.')
    bi_args.add_argument('--max_cont_length_cutoff', type=int, default=2250,
                         help='Maximum contig length. Caps NX filtering value.')
    bi_args.add_argument('--min_cont_length_cutoff_marker', type=int, default=2250,
                         help='Minimum contig length containing CheckM markers. Caps NX filtering value.')
    bi_args.add_argument('--max_cont_length_cutoff_marker', type=int, default=2250,
                         help='Maximum contig length containing CheckM markers. Caps NX filtering value.')
    bi_args.add_argument('--max_n_contigs', type=float, default=5.0e5,
                         help='Maximum number of contigs binny uses. If the number of available contigs after minimum '
                              'size filtering exceeds this, binny will increase the minimum size threshold until the '
                              'maximum is reached. Prevents use of excessive amounts of memory on large assemblies.'
                              'Default should ensure adequate performance, adjust according to available memory.')

    bi_args.add_argument('--max_marker_lineage_depth_lvl', type=int, default=2, choices=[0, 1, 2, 3, 4, 5, 6],
                         help='Maximum marker set lineage depth to check bin quality with:'
                              '0: \'domain\', 1: \'phylum\', 2: \'class\', 3: \'order\', 4: \'family\', 5: \'genus\', '
                              '6: \'species\'.')
    bi_args.add_argument('--distance_metric', type=str, default='manhattan',
                         choices=['cityblock', 'cosine', 'euclidean', 'haversine', 'l1', 'l2',
                                  'manhattan', 'nan_euclidean'],
                         help='Distance metric for opentSNE and HDBSCAN.')
    bi_args.add_argument('--max_iterations', type=int, default=50, help='Maximum number of binny iterations.')
    bi_args.add_argument('--hdbscan_epsilon_range', type=str, default='0.250,0.000',
                         help='Increasing the HDBSCAN cluster selection epsilon beyond 0.5 is not advised as it might'
                              ' massively increase run time, but it might help recover fragmented genomes that would be'
                              ' missed with lower settings.')
    bi_args.add_argument('--hdbscan_min_samples_range', type=str, default='1,5,10',
                         help='Adapted from the HDBSCAN manual: \'Measure of how conservative the clustering should be.'
                              ' With larger values, more points will be declared as noise, and clusters will be '
                              'restricted to progressively more dense areas.\'.')
    bi_args.add_argument('--include_depth_initial', action='store_true', default=False,
                         help='Use depth as additional dimension during the initial clustering.')
    bi_args.add_argument('--include_depth_main', action='store_true', default=False,
                         help='Use depth as additional dimension during the main clusterings.')
    bi_args.add_argument('--min_completeness', type=float, default=72.5,
                         help='Minimum value binny will lower completeness to while running.')
    bi_args.add_argument('--start_completeness', type=float, default=92.5,
                         help='Completeness threshold binny wilt begin with.')
    bi_args.add_argument('--purity', type=float, default=95,
                         help='Minimum purity for bins to be selected.')
    bi_args.add_argument('--write_contig_data', action='store_true', default=True,
                         help='Write all contig data to compressed tsv. Might create large file.')

    return parser


def binny_setup(args):
    # Try mamba, then conda
    if not run_cmd('which mamba').returncode:
        env_manager = 'mamba'
    elif not run_cmd('which conda').returncode:
        env_manager = 'mamba'
    else:
        raise Exception('Could not find Mamba or Conda in path.')

    # Make sure a compiler for cython is available
    if run_cmd('which gcc').returncode:
        print('Could not find gcc in path. Adding to conda environment.')
        run_cmd(f'{env_manager} install -c conda-forge gcc_linux-64 --yes')
        run_cmd('${CONDA_PREFIX}/etc/conda/activate.d/activate-binutils_linux-64.sh')
        run_cmd('${CONDA_PREFIX}/etc/conda/activate.d/activate-gcc_linux-64.sh')

    # Set up CheckM database
    if not os.path.exists(os.path.join(args.db_path, 'taxon_marker_sets_lineage_sorted.tsv')):
        print('Downloading CheckM data.')
        os.makedirs(args.db_path, exist_ok=True)
        checkm_data_url = 'https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz'
        download_file(checkm_data_url, output_folder=args.db_path + '/', stdout_file=None, retry_limit=10)

        checkm_tar = tarfile.open(os.path.join(args.db_path, 'checkm_data_2015_01_16.tar.gz'))
        checkm_tar.extract('./taxon_marker_sets.tsv', args.db_path)
        checkm_tar.extract('./pfam/tigrfam2pfam.tsv', args.db_path)
        checkm_tar.extract('./hmms/checkm.hmm', args.db_path)
        markers_df = pd.read_csv(os.path.join(args.db_path, 'taxon_marker_sets.tsv'), sep='\t', skipinitialspace=True,
                                 header=None)
        markers_df = markers_df.sort_values(markers_df.columns[2])
        markers_df.to_csv(os.path.join(args.db_path, 'taxon_marker_sets_lineage_sorted.tsv'), header=None, index=None,
                          sep='\t')
        print('Processing CheckM data.')
        prepCheckM.remove_unused_checkm_hmm_profiles(os.path.join(args.db_path, 'hmms/checkm.hmm'),
                                                     os.path.join(args.db_path, 'taxon_marker_sets.tsv'),
                                                     os.path.join(args.db_path, 'pfam/tigrfam2pfam.tsv'),
                                                     os.path.join(args.db_path, 'hmms'))
        if os.path.exists(os.path.join(args.db_path, 'checkm_data_2015_01_16.tar.gz')):
            os.remove(os.path.join(args.db_path, 'checkm_data_2015_01_16.tar.gz'))
        if os.path.exists(os.path.join(args.db_path, 'hmms/checkm.hmm')):
            os.remove(os.path.join(args.db_path, 'hmms/checkm.hmm'))
        if os.path.exists(os.path.join(args.db_path, 'taxon_marker_sets.tsv')) and os.path.exists(
                os.path.join(args.db_path, 'taxon_marker_sets_lineage_sorted.tsv')):
            os.remove(os.path.join(args.db_path, 'taxon_marker_sets.tsv'))

    # Prepare and install Mantis

    mantis_cfg = os.path.join(BINNY_DIR, 'config/binny_mantis.cfg')
    mantis_cfg_temp = os.path.join(BINNY_DIR, 'config/binny_mantis_template.cfg')
    prep_mantis_cfg = run_cmd(f'sed -e \'s|__PATH_TO_DB__|{args.db_path}|g\' {mantis_cfg_temp} > {mantis_cfg}')
    if prep_mantis_cfg.returncode:
        raise Exception(f'Failed to prepare Mantis config file.\n'
                        f'Error: {prep_mantis_cfg.stderr}')
    try:
        print('Checking Mantis setup.')
        check_installation(mantis_config=mantis_cfg, no_taxonomy=True)
    except InstallationCheckNotPassed:
        # Check for db setup success before general final check?
        setup_databases(mantis_config=mantis_cfg, no_taxonomy=True, chunk_size=1200, cores=args.threads)
        check_installation(mantis_config=mantis_cfg, no_taxonomy=True)


def run_binny_snakemake(args):
    # Check Snakefile
    if not os.path.exists(snakefile):
        msg = f'Error: cannot find Snakefile:\n{snakefile}\n'
        sys.stderr.write(msg)
        sys.exit(-1)

    # Check config file and setup Snakemake config.
    configs = [int_config]
    if args.config_file and not os.path.exists(args.config_file):
        msg = f'Error: cannot find config file:\n{args.config_file}\n'
        sys.stderr.write(msg)
        sys.exit(-1)
    elif args.config_file:
        configs.append(args.config_file)
        config = {}
    else:
        config = vars(args)
        if not args.tmp_dir:
            config['tmp_dir'] = f'{args.outputdir}/tmp'

    if args.use_cluster:
        clust_sub_cmd = args.cluster_submission_command
        clust_conf = args.scheduler_preset_path
    else:
        clust_sub_cmd = None
        clust_conf = None

    status = snakemake.snakemake(snakefile, configfiles=configs, cores=args.threads, config=config,
                                 printshellcmds=True, dryrun=args.dry_run, forceall=args.force, cluster=clust_sub_cmd,
                                 cluster_config=clust_conf, jobname=args.job_name, unlock=args.unlock)

    if args.report and not args.dry_run:
        report_file = os.path.join(args.outputdir, 'report.html')
        snakemake.snakemake(snakefile, configfiles=configs, config=config, report=report_file)

    if status:  # translate "success" into shell exit code of 0
        return 0
    return 1

def main():
    parser = arg_parser()
    arguments = parser.parse_args()
    binny_setup(arguments)
    run_binny_snakemake(arguments)

if __name__ == '__main__':
    main()