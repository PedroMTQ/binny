#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 22 10:50:35 2020

@author: oskar.hickl
"""

mg_depth_file = snakemake.input['mgdepth']
assembly = snakemake.input['assembly']
annot_file = snakemake.input['gff']
functions = snakemake.params['py_functions']
min_purity = snakemake.params['purity']
min_completeness = snakemake.params['completeness']
kmers = snakemake.params['kmers']
min_contig_length = snakemake.params['cutoff']
threads = snakemake.threads
log = snakemake.log[0]

import sys
sys.path.append(functions)
from binny_functions import *

sys.stdout = open(log, 'w')

# Load assembly and mask rRNAs and CRISPR arrays
contig_list = [[key] + [val] for key, val in load_fasta(assembly).items() if len(val) >= int(min_contig_length)]
contig_rrna_crispr_region_dict = gff2low_comp_feature_dict(annot_file)
mask_rep_featrues(contig_rrna_crispr_region_dict, contig_list)

# Get length normalized k-mer frequencies.
kmer_sizes = kmers.split(',')
start = timer()
kfreq_array = get_contig_kmer_matrix2(contig_list, kmer_sizes, threads)
end = timer()
print('K-mer frequency matrix created in {0}s.'.format(end - start))

# Standardize k-mer freq data
X = np.array([c_kfreq[1:] for c_kfreq in kfreq_array[1:]])
X_contigs = [c_kfreq[0] for c_kfreq in kfreq_array[1:]]
scaler = StandardScaler().fit(X)
X_scaled = scaler.transform(X)

# Manifold learning and dimension reduction.
print('Running manifold learning and dimension reduction')
n_dim = 3
perp = get_perp(len(X_contigs))
# print('Perplexity: {0}.'.format(perp))
X_embedded = umap.UMAP(n_neighbors=perp, min_dist=0.0, n_components=n_dim, random_state=0, densmap=False, verbose=True).fit_transform(X_scaled)

# Create 2/3D coordinate df.
if n_dim == 3:
    coord_df = pd.DataFrame(data=X_embedded, index=None, columns=['x', 'y', 'z'])
    coord_df['contig'] = X_contigs
    coord_df = coord_df[['contig', 'x', 'y', 'z']]
    coord_df.to_csv('intermediary/contig_coordinates.tsv', sep='\t', index=False, header=False)
    coords_file = 'intermediary/contig_coordinates.tsv'
elif n_dim == 2:
    coord_df = pd.DataFrame(data=X_embedded, index=None, columns=['x', 'y'])
    coord_df['contig'] = X_contigs
    coord_df = coord_df[['contig', 'x', 'y']]
    coord_df.to_csv('intermediary/contig_coordinates.tsv', sep='\t', index=False, header=False)
    coords_file = 'intermediary/contig_coordinates.tsv'

# Load data
contig_data_df = load_and_merge_cont_data(annot_file, mg_depth_file, coords_file, sne_dims=n_dim)
# Write contig data to file
contig_data_df.to_csv('intermediary/contig_data.tsv', sep='\t', index=False)

# Find bins
good_bins, final_init_clust_dict = binny_iterate(contig_data_df, threads, min_purity, min_completeness)

# Write final table
final_clust_dict = {**good_bins, **final_init_clust_dict}
final_clust_df = cluster_df_from_dict(final_clust_dict)
final_clust_df.to_csv('final_contigs2clusters.tsv', sep='\t', index=False)

final_clust_contig_df = contig_df_from_cluster_dict(final_clust_dict)

# Plot final clustering
conditions = [(final_clust_contig_df['purity'] >= min_purity / 100) & (
            final_clust_contig_df['completeness'] >= min_completeness / 100),
              (final_clust_contig_df['purity'] < min_purity / 100) | (
                          final_clust_contig_df['completeness'] < min_completeness / 100)]
values = [final_clust_contig_df['cluster'], 'N']
final_clust_contig_df['above_thresh'] = np.select(conditions, values)
final_clust_contig_df = contig_data_df.merge(final_clust_contig_df, how='outer', on='contig',
                                                 suffixes=(None, '_y'))
final_clust_contig_df['above_thresh'] = final_clust_contig_df['above_thresh'].fillna('N')

write_scatterplot(final_clust_contig_df, 'final_scatter_plot.pdf',
                  final_clust_contig_df['above_thresh'])

# Check if something went wrong and contigs were duplicated.
all_contigs = []
for bin in good_bins:
    all_contigs.extend(good_bins[bin]['contigs'])
if len(all_contigs) != len(set(all_contigs)):
    print('WARNING: Duplicate contigs in bins found!')

# Write bin fastas.
write_bins(good_bins, assembly, min_comp=int(min_completeness), min_pur=int(min_purity), bin_dir='bins')

print('Run finished.')
sys.stdout.close()
