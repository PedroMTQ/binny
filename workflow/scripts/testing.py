binny_functions_path = '/Users/oskar.hickl/local_tools/binny_fork/binny/workflow/scripts/'

import sys
sys.path.append(binny_functions_path)
from binny_functions import *
import khmer

import itertools
import re
from timeit import default_timer as timer

########################################################################################################################
import os
os.chdir('/Users/oskar.hickl/binny_bench/binner_data/Binny/binny_outputs/2017.12.04_18.45.54_sample_5')
########################################################################################################################
# cluster_dir = snakemake.input['outdir']
# binny_dir = snakemake.params['binnydir']
binnydir = 'intermediary/'
# mg_depth_file = snakemake.input['mgdepth']
mg_depth_file = 'intermediary/assembly.contig_depth.txt'
# assembly = snakemake.input['assembly']
assembly = 'assembly.fa'
# assembly_cut = 'intermediary/assembly.cut.fa'
# coords_file = snakemake.input['vizbin']
# coords_file = 'vizbin.with-contig-names.points'
# annot_file = snakemake.input['gff']
annot_file = 'intermediary/annotation_CDS_RNA_hmms.gff'
# gff_filt = 'intermediary/annotation.filt.gff'
# function_script =snakemake.params['plot_functions']
# pk = snakemake.config['binning']['binny']['pk']
# pk = 10
# all_outputs = snakemake.output
# gff_db = 'intermediary/annotation_CDS_RNA_hmms.gff.db'
# gff_db = 'intermediary/annotation_CDS_RNA_hmms.sqlt'

min_completeness = 40
min_purity = 90
threads = 10


# # Load assembly and mask rRNAs and CRISPR arrays
# contig_list = [[key] + [val] for key, val in load_fasta(assembly).items() if len(val) >= 500]
# contig_rrna_crispr_region_dict = gff2low_comp_feature_dict(gff_filt)
# mask_rep_featrues(contig_rrna_crispr_region_dict, contig_list)
#
# # Get length normalized k-mer frequencies.
# kmer_sizes = list(range(2, 5))
# start = timer()
# kfreq_array = get_contig_kmer_matrix2(contig_list, kmer_sizes, threads)
# end = timer()
# print(end - start)
#
# # Standardize k-mer freq data
# X = np.array([c_kfreq[1:] for c_kfreq in kfreq_array[1:]])
# X_contigs = [c_kfreq[0] for c_kfreq in kfreq_array[1:]]
# scaler = StandardScaler().fit(X)
# X_scaled = scaler.transform(X)

# n_comp = None
# sum_var_exp = 0
# n_tries = 0
# if len(X_scaled[0]-1) < 25:
#     n_comp = len(X_scaled[0]-1)
# else:
#     n_comp = 25
#
# pca = PCA(n_components=n_comp)
# transformer = pca.fit(X_scaled)
# sum_var_exp = sum(pca.explained_variance_ratio_)
# print(n_comp, sum(pca.explained_variance_ratio_), pca.explained_variance_ratio_)
#
# while sum_var_exp <= 0.90 and n_tries <= 100:
#     pca = PCA(n_components=n_comp)
#     transformer = pca.fit(X_scaled)
#     sum_var_exp = sum(pca.explained_variance_ratio_)
#     # print(pca.explained_variance_ratio_)
#     # print(pca.singular_values_)
#     n_comp += 25
#     n_tries += 1
# print(n_comp, sum(pca.explained_variance_ratio_), pca.explained_variance_ratio_)
# X_pca = transformer.transform(X_scaled)

# n_dim = 3
# perp = get_perp(len(X_contigs))
# print('Perplexity: {0}.'.format(perp))
#
# X_embedded = umap.UMAP(n_neighbors=perp, min_dist=0.0, n_components=n_dim, random_state=0, densmap=False, verbose=True).fit_transform(X_scaled)

# learning_rate = 10  # 2
# early_exaggeration = int(learning_rate / 20) + 1  # 1
# X_embedded = TSNE(n_components=n_dim, perplexity=perp, learning_rate=learning_rate, early_exaggeration=early_exaggeration,
#                   metric='euclidean', init='pca', n_jobs=threads, verbose=2, random_state=0).fit_transform(X_pca)

# sns.set_style("whitegrid", {'axes.grid' : False})
#
# fig = plt.figure(figsize=(7, 6))
#
#
# if n_dim == 3:
#     df = pd.DataFrame(X_embedded, columns=['x', 'y', 'z'])
#     df['contig'] = X_contigs
#     # ax = fig.add_subplot(111, projection='3d')
#     ax = Axes3D(fig)
#     if n_comp:
#         ax.set_title("N_DIMS: {0}, PERP: {1}, LR: {2}, EE: {3}, PCA_COMP: {4}".format(n_dim, perp, learning_rate, early_exaggeration, n_comp))
#     else:
#         ax.set_title("N_DIMS: {0}, PERP: {1}, LR: {2}, EE: {3}".format(n_dim, perp, learning_rate, early_exaggeration))
#     ax.patch.set_facecolor('orange')
#     x = df['x']
#     y = df['y']
#     z = df['z']
#     ax.scatter(x, y, z, s=0.05, alpha=1)
# elif n_dim == 2:
#     df = pd.DataFrame(X_embedded, columns=['x', 'y'])
#     df['contig'] = X_contigs
#     if n_comp:
#         sns.scatterplot(data=df, x="x", y="y", sizes=0.1, s=1).set_title("N_DIMS: {0}, PERP: {1}, LR: {2}, EE: {3}, PCA_COMP: {4}".format(n_dim,
#                                                                                                                            perp,
#                                                                                                                            learning_rate,
#                                                                                                                            early_exaggeration,
#                                                                                                                            n_comp))
#     else:
#         sns.scatterplot(data=df, x="x", y="y", sizes=0.1, s=1).set_title(
#             "N_DIMS: {0}, PERP: {1}, LR: {2}, EE: {3}".format(n_dim,
#                                                                              perp,
#                                                                              learning_rate,
#                                                                              early_exaggeration))
# plt.show()

########################################################################################################################

# Create 2/3D coordinate df.
# if n_dim == 3:
#     coord_df = pd.DataFrame(data=X_embedded, index=None, columns=['x', 'y', 'z'])
#     coord_df['contig'] = X_contigs
#     coord_df = coord_df[['contig', 'x', 'y', 'z']]
#     coord_df.to_csv('contig_coordiantes_py_sne.tsv', sep='\t', index=False, header=False)
#     coords_file = 'contig_coordiantes_py_sne.tsv'
# elif n_dim == 2:
#     coord_df = pd.DataFrame(data=X_embedded, index=None, columns=['x', 'y'])
#     coord_df['contig'] = X_contigs
#     coord_df = coord_df[['contig', 'x', 'y']]
#     coord_df.to_csv('contig_coordiantes_py_sne.tsv', sep='\t', index=False, header=False)
#     coords_file = 'contig_coordiantes_py_sne.tsv'
########################################################################################################################

# # Load data
# contig_data_df = load_and_merge_cont_data(annot_file, mg_depth_file, coords_file, sne_dims=n_dim)
# # Write contig data to file
# contig_data_df.to_csv('contig_data_v13.tsv', sep='\t', index=False)
#
# # Find bins
# good_bins = binny_iterate(contig_data_df, threads, min_purity, min_completeness)
#
# # Write final table
#
# final_clust_dict = {**good_clusters, **init_clust_dict}
# final_clust_df = cluster_df_from_dict(final_clust_dict)
# final_clust_df.to_csv('final_{0}_contigs2clusters_{0}.tsv'.format(n_iterations), sep='\t', index=False)
#
# final_clust_contig_df = contig_df_from_cluster_dict(final_clust_dict)
#
# # Plot final clustering
# conditions = [(final_clust_contig_df['purity'] >= min_purity / 100) & (
#             final_clust_contig_df['completeness'] >= min_completeness / 100),
#               (final_clust_contig_df['purity'] < min_purity / 100) | (
#                           final_clust_contig_df['completeness'] < min_completeness / 100)]
# values = [final_clust_contig_df['cluster'], 'N']
# final_clust_contig_df['above_thresh'] = np.select(conditions, values)
# final_clust_contig_df = contig_data_df.merge(final_clust_contig_df, how='outer', on='contig',
#                                                  suffixes=(None, '_y'))
# final_clust_contig_df['above_thresh'] = final_clust_contig_df['above_thresh'].fillna('N')
#
# write_scatterplot(final_clust_contig_df, 'final_{0}_final_scatter_plot_v13_{0}.pdf'.format(n_iterations),
#                   final_clust_contig_df['above_thresh'])
#
# # Check if something went wrong and contigs were duplicated.
# all_contigs = []
# for bin in good_clusters:
#     all_contigs.extend(good_clusters[bin]['contigs'])
# if len(all_contigs) != len(set(all_contigs)):
#     print('WARNING: Duplicate contigs found!')
#
# # Write bin fastas.
# write_bins(good_clusters, assembly, min_comp=int(min_completeness), min_pur=int(min_purity), bin_dir='bins_v13')


########################################################################################################################
### Write original optimized binny out
########################################################################################################################
# path2tsv = '/Users/oskar.hickl/binny_bench/binner_data/Binny/binny_outputs/2017.12.04_18.45.54_sample_5_50/contigs2bin.tsv'
# org_binny_df = pd.read_csv(path2tsv, sep='\t')
#
# org_contig_data_df = load_and_merge_cont_data(annot_file, mg_depth_file, 'vizbin.with-contig-names.points', sne_dims=n_dim)
# data_df_2d = org_binny_df.merge(org_contig_data_df[['contig', 'x', 'y']], how='left', on='contig')
#
# write_scatterplot(data_df_2d, 'original_binny_out.pdf', data_df_2d['cluster'])
########################################################################################################################


# Load assembly and mask rRNAs and CRISPR arrays
contig_list = [[key] + [val] for key, val in load_fasta(assembly).items() if len(val) >= 500]
contig_rrna_crispr_region_dict = gff2low_comp_feature_dict(annot_file)
mask_rep_featrues(contig_rrna_crispr_region_dict, contig_list)

# Get length normalized k-mer frequencies.
kmer_sizes = list(range(2, 5))
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
final_clust_df.to_csv('final_contigs2clusters_sc.tsv', sep='\t', index=False)

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

write_scatterplot(final_clust_contig_df, 'final_scatter_plot_sc.pdf',
                  final_clust_contig_df['above_thresh'])

# Check if something went wrong and contigs were duplicated.
all_contigs = []
for bin in good_bins:
    all_contigs.extend(good_bins[bin]['contigs'])
if len(all_contigs) != len(set(all_contigs)):
    print('WARNING: Duplicate contigs in bins found!')

# Write bin fastas.
write_bins(good_bins, assembly, min_comp=int(min_completeness), min_pur=int(min_purity), bin_dir='bins_sc')

print('Run finished.')
