"""
Created on Wed Feb 22 10:50:35 2020

@author: oskar.hickl
"""

from pathlib import Path
from joblib import parallel_backend, Parallel, delayed
import numpy as np
import pandas as pd
import seaborn as sns
from sklearn.neighbors import KNeighborsClassifier
from sklearn.cluster import DBSCAN
from sklearn.cluster import OPTICS
from unidip import UniDip
from timeit import default_timer as timer

import khmer
import itertools
import re
from timeit import default_timer as timer


from sklearn.manifold import TSNE
from sklearn.preprocessing import RobustScaler, StandardScaler
from sklearn.decomposition import PCA
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

import umap
import hdbscan

def unify_multi_model_genes(gene):
    # Dict to unify genes with multiple models
    hmm_dict = {'TIGR00388': 'glyS', 'TIGR00389': 'glyS', 'TIGR00471': 'pheT', 'TIGR00472': 'pheT', 'TIGR00408': 'proS',
                'TIGR00409': 'proS', 'TIGR02386': 'rpoC', 'TIGR02387': 'rpoC'}
    # Replace gene with unified name if it has multiple models
    if gene in hmm_dict:
        gene = hmm_dict.get(gene)
    return gene


def gff2ess_gene_df(annotation_file):
    contig_dict = {}
    # Get all info for each line
    with open(annotation_file, 'r') as af:
        for line in af:
            line = line.strip(' \t\n').split('\t')
            contig, contig_start, contig_end, attributes = line[0], line[3], line[4], line[8].split(';')
            # Check if feature has no start and/or stop position or length = 0
            if not contig_start or not contig_end or int(contig_start)-int(contig_end) == 0:
                continue
            # Only add line if feature was annotated as essential gene
            for attribute in attributes:
                attribute_name, attribute_id = attribute.split('=')[0], attribute.split('=')[1]
                if attribute_name == 'essential' and attribute_id:
                    if not contig_dict.get(contig):
                        contig_dict[contig] = [unify_multi_model_genes(attribute_id)]
                    else:
                        contig_dict[contig].append(unify_multi_model_genes(attribute_id))
    # Create list of lists to join dict items as string and build data frame
    annot_cont_line_list = [[contig, ','.join(contig_dict.get(contig))] for contig in contig_dict]
    annot_cont_df = pd.DataFrame(annot_cont_line_list, columns=['contig', 'essential'])
    return annot_cont_df


def load_and_merge_cont_data(annot_file, depth_file, sne_coords, sne_dims=2):
    print('Reading Prokka annotation.')
    annot_df = gff2ess_gene_df(annot_file)
    print('Reading BH-tSNE coordinates.')
    if sne_dims == 3:
        coord_df = pd.read_csv(sne_coords, sep='\t', names=['contig', 'x', 'y', 'z'])
    else:
        coord_df = pd.read_csv(sne_coords, sep='\t', names=['contig', 'x', 'y'])
    print('Reading average depth.')
    depth_df = pd.read_csv(depth_file, sep='\t', names=['contig', 'depth'])
    print('Merging data')
    # Merge annotation data and VizBin coords first, keeping only contigs with coords, then merge with depth data
    cont_data_df = annot_df.merge(coord_df, how='right', on='contig').merge(depth_df, how='inner', on='contig').sort_values(by='contig')
    return cont_data_df


def knn_sne_coords(contig_info_df, knn_sne_pk, sne_dims=2):
    if knn_sne_pk < 2 * 2:
        knn_sne_pk = 2 * 2
    if knn_sne_pk >= contig_info_df['contig'].size:
        knn_sne_pk = 2
    knn_c = KNeighborsClassifier(n_neighbors=knn_sne_pk, weights='distance')
    if sne_dims == 3:
        coords_df = contig_info_df.loc[:, ['x', 'y', 'z']]
    else:
        coords_df = contig_info_df.loc[:, ['x', 'y']]
    knn_c.fit(coords_df, contig_info_df.loc[:, 'contig'])
    # sort distance of neighbouring points
    try:
        skn = pd.Series(sorted([row[-1] for row in knn_c.kneighbors()[0]]))
    except ValueError:
        print(knn_sne_pk, contig_info_df, coords_df, contig_info_df.loc[:, 'contig'])
        raise Exception
    # calculate running standard deviation between 10 neighbouring points
    sdkn = skn.rolling(10, center=True, min_periods=0).std()
    # find the first jump in distances at the higher end
    try:
        est = sorted([e for i, e in zip(sdkn, skn) if i > sdkn.quantile(.975) and e > skn.median()])[0]
    except IndexError:
        est = None
    return est


def contig_df2cluster_dict(contig_info_df, dbscan_labels, use_noise=False):
    # This function is garbage. From pandas df to dict to df. Need to improve asap
    if use_noise:
        use_noise = 'noise_stays'
    else:
        use_noise = '-1'
    cluster_dict = {}
    dims = [dim for dim in ['x', 'y', 'z'] if dim in contig_info_df.columns]
    # print(dims)
    contig_info_df = contig_info_df.loc[:, ['contig', 'essential', 'depth'] + dims]
    contig_info_df['cluster'] = dbscan_labels
    tmp_contig_dict = contig_info_df.fillna('non_essential').set_index('contig').to_dict('index')
    for contig in tmp_contig_dict:
        contig_cluster = shorten_cluster_names(str(tmp_contig_dict.get(contig, {}).get('cluster')))
        contig_essential = tmp_contig_dict.get(contig, {}).get('essential')
        contig_depth = tmp_contig_dict.get(contig, {}).get('depth')
        contig_x = tmp_contig_dict.get(contig, {}).get('x')
        contig_y = tmp_contig_dict.get(contig, {}).get('y')
        if tmp_contig_dict.get(contig, {}).get('z'):
            contig_z = tmp_contig_dict.get(contig, {}).get('z')
        if not cluster_dict.get(contig_cluster) and use_noise not in contig_cluster:
            cluster_dict[contig_cluster] = {'essential': [contig_essential],  # .split(','),
                                            'depth': [contig_depth],
                                            'contigs': [contig],
                                            'x': [contig_x],
                                            'y': [contig_y]}
            if tmp_contig_dict.get(contig, {}).get('z'):
                cluster_dict[contig_cluster]['z'] = [contig_z]
        elif use_noise not in contig_cluster:
            if contig_essential:
                cluster_dict.get(contig_cluster, {}).get('essential').append(contig_essential)  # .extend(contig_essential.split(','))
            cluster_dict.get(contig_cluster, {}).get('depth').append(contig_depth)
            cluster_dict.get(contig_cluster, {}).get('contigs').append(contig)
            cluster_dict.get(contig_cluster, {}).get('x').append(contig_x)
            cluster_dict.get(contig_cluster, {}).get('y').append(contig_y)
            if tmp_contig_dict.get(contig, {}).get('z'):
                cluster_dict.get(contig_cluster, {}).get('z').append(contig_z)
    return cluster_dict


def dbscan_cluster(contig_data_df, pk=None, include_depth=False, n_jobs=1):
    dims = [dim for dim in ['x', 'y', 'z'] if dim in contig_data_df.columns]
    if not include_depth:
        dim_df = contig_data_df.loc[:, dims].to_numpy(dtype=np.float64)
    else:
        dim_df = contig_data_df.loc[:, dims + ['depth']].to_numpy(dtype=np.float64)
    if not pk:
        pk = int(np.log(contig_data_df['contig'].size))
        if pk < len(dim_df) * 2:
            pk = len(dim_df) * 2
    # Get reachability distance estimate
    est = knn_sne_coords(contig_data_df, pk)
    # Run parallelized dbscan
    # print('Running dbscan.')
    with parallel_backend('threading'):
        dbsc = DBSCAN(eps=est, min_samples=pk, n_jobs=n_jobs).fit(dim_df)
    cluster_labels = dbsc.labels_
    if len(set(cluster_labels)) == 1 and str(list(set(cluster_labels))[0]) == '-1':
        cluster_dict = contig_df2cluster_dict(contig_data_df, cluster_labels, use_noise=True)
    else:
        cluster_dict = contig_df2cluster_dict(contig_data_df, cluster_labels)
    return cluster_dict, cluster_labels


def optics_cluster(contig_data_df, min_samples=None, include_depth=False, n_jobs=1):
    dims = [dim for dim in ['x', 'y', 'z'] if dim in contig_data_df.keys()]
    if not include_depth:
        dim_df = contig_data_df.loc[:, dims].to_numpy(dtype=np.float64)
    else:
        dim_df = contig_data_df.loc[:, dims + ['depth']].to_numpy(dtype=np.float64)
    if not min_samples:
        min_samples = int(np.log(contig_data_df['contig'].size))
    # print('Running OPTICS.')
    clustering = OPTICS(min_samples=min_samples, n_jobs=n_jobs).fit(dim_df)
    cluster_labels = clustering.labels_
    if len(set(cluster_labels)) == 1 and str(list(set(cluster_labels))[0]) == '-1':
        cluster_dict = contig_df2cluster_dict(contig_data_df, cluster_labels, use_noise=True)
    else:
        cluster_dict = contig_df2cluster_dict(contig_data_df, cluster_labels)
    return cluster_dict, cluster_labels


def hdbscan_cluster(contig_data_df, pk=None, include_depth=False, n_jobs=1):
    dims = [dim for dim in ['x', 'y', 'z'] if dim in contig_data_df.columns]
    if not include_depth:
        dim_df = contig_data_df.loc[:, dims].to_numpy(dtype=np.float64)
    else:
        dim_df = contig_data_df.loc[:, dims + ['depth']].to_numpy(dtype=np.float64)
    if not pk:
        pk = get_perp(contig_data_df['contig'].size)
        # pk = int(np.log(contig_data_df['contig'].size))
        if pk < len(dims) * 2:
            pk = len(dims) * 2
    with parallel_backend('threading'):
        hdbsc = hdbscan.HDBSCAN(core_dist_n_jobs=n_jobs, min_cluster_size=pk).fit(dim_df)  # min_cluster_size=15
    cluster_labels = hdbsc.labels_

    while len(set(cluster_labels)) == 1 and str(list(set(cluster_labels))[0]) == '-1' and pk >= len(dims) * 2:
        pk = int(pk * 0.5)
        print('HDBSCAN found only noise trying with lower min_cluster_size={0}.'.format(pk))
        with parallel_backend('threading'):
            hdbsc = hdbscan.HDBSCAN(core_dist_n_jobs=n_jobs, min_cluster_size=pk).fit(dim_df)  # min_cluster_size=15
        cluster_labels = hdbsc.labels_

    if len(set(cluster_labels)) == 1 and str(list(set(cluster_labels))[0]) == '-1':
        cluster_dict = contig_df2cluster_dict(contig_data_df, cluster_labels, use_noise=True)
    else:
        cluster_dict = contig_df2cluster_dict(contig_data_df, cluster_labels)
    return cluster_dict, cluster_labels


def run_initial_scan(contig_data_df, initial_cluster_mode, dbscan_threads, pk=None, include_depth=False):
    if initial_cluster_mode == 'DBSCAN' or not initial_cluster_mode:
        if not pk:
            pk = get_perp(contig_data_df['contig'].size)
            # pk = int(np.log(contig_data_df['contig'].size))
        print('Running initial scan with DBSCAN and min samples of: {0}'.format(str(pk)))
        # Run parallelized dbscan
        first_clust_dict, labels = dbscan_cluster(contig_data_df, pk=pk, n_jobs=dbscan_threads,
                                                  include_depth=include_depth)
    elif initial_cluster_mode == 'OPTICS':
        if not pk:
            # pk = int(pow(np.log(contig_data_df['contig'].size), 2) / 2)
            # pk = int(np.log(contig_data_df['contig'].size))
            # pk = int(np.log((contig_data_df['contig'].size) ** (np.log(contig_data_df['contig'].size) / np.log(10))))
            pk = get_perp(contig_data_df['contig'].size)
        print('Running initial scan with OPTICS and min samples of: {0}'.format(str(pk)))
        # Run OPTICS
        first_clust_dict, labels = optics_cluster(contig_data_df, min_samples=pk, n_jobs=dbscan_threads,
                                              include_depth=include_depth)
    elif initial_cluster_mode == 'HDBSCAN':
        if not pk:
            pk = get_perp(contig_data_df['contig'].size)
            # pk = int(np.log(contig_data_df['contig'].size))
        print('Running initial scan with HDBSCAN and min samples of: {0}'.format(str(pk)))
        # Run parallelized dbscan
        first_clust_dict, labels = hdbscan_cluster(contig_data_df, pk=pk, n_jobs=dbscan_threads,
                                                  include_depth=include_depth)

    return first_clust_dict, labels


def cluster_df_from_dict(cluster_dict):
    cluster_df = pd.DataFrame()
    if len(set(list(cluster_dict.keys()))) == 1 and str(list(set(list(cluster_dict.keys())))[0]) == '-1':
        cluster_df['cluster'] = [cluster for cluster in cluster_dict]
    else:
        cluster_df['cluster'] = [cluster for cluster in cluster_dict if not cluster == '-1']
    for metric in ['contigs', 'essential']:
        metric_list = []
        metric_uniq_list = None
        if metric == 'essential':
            metric_uniq_list = []
        for cluster in cluster_dict:
            essential_list = [i for e in cluster_dict.get(cluster, {}).get(metric) for i in e.split(',')
                              if not e == 'non_essential']
            metric_list.append(len(essential_list))
            if metric == 'essential':
                metric_uniq_list.append(len(set(essential_list)))
        cluster_df[metric] = metric_list
        if metric_uniq_list:
            try:
                cluster_df['unique_' + metric] = metric_uniq_list
            except ValueError:
                print(metric_uniq_list, metric, cluster_df, cluster_dict)
                raise Exception
    cluster_df['completeness'] = round(cluster_df['unique_essential'] / 115, 2)
    cluster_df['purity'] = round(cluster_df['unique_essential'] / cluster_df['essential'], 2)
    return cluster_df.fillna(0)


def cluster_dict2bin_dict(cluster_dict, cluster_df):
    cluster_df = cluster_df.set_index('cluster')
    bin_dict = {}
    for cluster in cluster_dict:
        cluster_purity = str(int(cluster_df.loc[cluster, 'purity'] * 100))
        cluster_completeness = str(int(cluster_df.loc[cluster, 'completeness'] * 100))
        bin_name = '_'.join(['binny'] + [cluster] + [cluster_purity] + [cluster_completeness])
        bin_dict[bin_name] = cluster_dict.get(cluster, {}).get('contigs')
    return bin_dict


def sort_cluster_dict_data_by_depth(cluster_dict):
    sorted_dict = {}
    dims = [dim for dim in ['x', 'y', 'z'] if dim in list(cluster_dict[list(cluster_dict.keys())[0]].keys())]
    for cluster in cluster_dict:
        sorted_dict[cluster] = {}
        deps = np.array(cluster_dict[cluster]['depth']).argsort()
        for metric in ['depth', 'contigs', 'essential'] + dims:
            metric_np = np.array(cluster_dict[cluster][metric])
            sorted_dict[cluster][metric] = metric_np[deps]
    return sorted_dict


def gather_cluster_data(cluster, cluster_dict):
    cluster_essential_genes = [gene for genes in cluster_dict.get(cluster, {}).get('essential')
                               for gene in genes.split(',') if not gene == 'non_essential' ]
    cluster_unique_essential_genes = set(cluster_essential_genes)
    if cluster_essential_genes:
        cluster_purity = round(len(cluster_unique_essential_genes) / len(cluster_essential_genes), 2)
        cluster_completeness = round(len(cluster_unique_essential_genes) / 115, 2)
    else:
        cluster_purity = 0
        cluster_completeness = 0
    cluster_data = [cluster_dict.get(cluster, {}).get('contigs'), cluster_dict.get(cluster, {}).get('essential'),
                    cluster_dict.get(cluster, {}).get('depth'), cluster_dict.get(cluster, {}).get('x'),
                    cluster_dict.get(cluster, {}).get('y'), cluster_purity, cluster_completeness]
    if 'z' in list(cluster_dict[cluster].keys()):
        cluster_data = cluster_data[:4] + [cluster_dict.get(cluster, {}).get('z')] + cluster_data[4:]
    return cluster_data


def create_new_clusters(cluster, intervals, clust_dat):
    new_clusters = {}
    assigned_indices = []
    bin_id = 1
    sum_sub_cluster_indices = 0
    for interval in intervals:
        cluster_id = cluster + '.' + str(bin_id)
        new_clusters[cluster_id] = {'depth': clust_dat[2][interval[0]:interval[1]],
                                    'essential': clust_dat[1][interval[0]:interval[1]],
                                    'contigs': clust_dat[0][interval[0]:interval[1]],
                                    'x': clust_dat[3][interval[0]:interval[1]],
                                    'y': clust_dat[4][interval[0]:interval[1]]}
        if len(clust_dat) == 8:
            new_clusters[cluster_id]['z'] = clust_dat[5][interval[0]:interval[1]]
        bin_id += 1
        # remove indices from leftover list
        assigned_indices.extend(range(interval[0], interval[1]))
        sum_sub_cluster_indices += len(clust_dat[2][interval[0]:interval[1]])
    # Write leftover 'cluster'
    cluster_id = cluster + '.L'
    new_clusters[cluster_id] = {'depth': np.array([clust_dat[2][index] for index, i in enumerate(clust_dat[2].tolist())
                                                   if index not in assigned_indices]),
                                'essential': np.array(
                                    [clust_dat[1][index] for index, i in enumerate(clust_dat[1].tolist())
                                     if index not in assigned_indices]),
                                'contigs': np.array(
                                    [clust_dat[0][index] for index, i in enumerate(clust_dat[0].tolist())
                                     if index not in assigned_indices]),
                                'x': np.array([clust_dat[3][index] for index, i in enumerate(clust_dat[3].tolist())
                                               if index not in assigned_indices]),
                                'y': np.array([clust_dat[4][index] for index, i in enumerate(clust_dat[4].tolist())
                                               if index not in assigned_indices])}
    if len(clust_dat) == 8:
        new_clusters[cluster_id]['z'] = np.array([clust_dat[5][index] for index, i in enumerate(clust_dat[5].tolist())
                                               if index not in assigned_indices])
    bin_id += 1
    return new_clusters


# def get_sub_clusters(cluster, cluster_dict, threads_for_dbscan, purity_threshold=0.9, completeness_threshold=0.4, alpha=0.05, pk=None, cluster_mode=None, include_depth=False):
#     # Create dict with just cluster and sort again, to ensure order by depth is
#     cluster_dict = sort_cluster_dict_data_by_depth({cluster: cluster_dict[cluster]})
#     # All data needed stored in list with following order:
#     # 0:contigs, 1:genes, 2:depth, 3:x, 4:y, 5:purity, 6:completeness
#     clust_dat = gather_cluster_data(cluster, cluster_dict)
#     if clust_dat[5] < purity_threshold and isinstance(clust_dat[5], float) and clust_dat[6] > completeness_threshold:
#         print('Cluster {0} below purity threshold of {1} with {2} and above completeness threshold of {3} with {4}. '
#               'Attempting to split.'.format(shorten_cluster_names(cluster), purity_threshold, clust_dat[5], completeness_threshold, clust_dat[6]))
#
#         # intervals = [1]
#         # start = timer()
#         # while len(intervals) == 1 and alpha <= 0.1:
#         #     intervals = UniDip(clust_dat[2], alpha=alpha, ntrials=1000, mrg_dst=1).run()
#         #     alpha += 0.025
#         #     if len(intervals) > 1:
#         #         end = timer()
#         #         # print('Found {0} depth sub-clusters in cluster {1} with alpha of {2} after {3}s.'.format(
#         #         #     len(intervals), shorten_cluster_names(cluster), round(alpha, 3), int(end - start)))
#         #         new_clusters = create_new_clusters(cluster, intervals, clust_dat)
#         #         return [new_clusters, cluster]
#         # end = timer()
#         # print('Failed to find depth sub-bclusters for {0} after {1}s.'.format(shorten_cluster_names(cluster), int(end - start)))
#         # initialise some stuff
#
#         cluster_contig_df = contig_df_from_cluster_dict(cluster_dict)
#         min_dims = 2 * 2
#         if min_dims > cluster_contig_df['contig'].size:
#             min_dims = 2
#         if not pk:
#             pk = int(np.log(cluster_contig_df['contig'].size))
#         if pk < min_dims:
#             pk = min_dims
#         new_clusters_labels = [1]
#         if cluster_mode == 'DBSCAN' or not cluster_mode:
#             print('Trying with DBSCAN.')
#             while cluster_contig_df['contig'].size < pk and pk >= min_dims:
#                 pk = int(pk * 0.75)
#             dbscan_tries = 0
#             print('Working on {0}.'.format(shorten_cluster_names(cluster)))
#             while len(set(new_clusters_labels)) == 1 and min_dims <= pk < cluster_contig_df['contig'].size:
#                 dbscan_tries += 1
#                 cluster_est = knn_sne_coords(cluster_contig_df, pk)
#                 while not cluster_est and pk >= min_dims:
#                     pk = int(pk * 0.75)
#                     if pk >= min_dims:
#                         cluster_est = knn_sne_coords(cluster_contig_df, pk)
#                 if cluster_est and cluster_contig_df['contig'].size > pk >= min_dims:
#                     df_vizbin_coords = cluster_contig_df.loc[:, ['x', 'y']].to_numpy(dtype=np.float64)
#                     with parallel_backend('threading'):
#                         dbsc = DBSCAN(eps=cluster_est, min_samples=pk, n_jobs=threads_for_dbscan).fit(df_vizbin_coords)
#                     new_clusters_labels = dbsc.labels_
#                     if len(set(new_clusters_labels)) > 1:
#                         new_cluster_names = {item: cluster + '.' + str(index + 1) for index, item in
#                                              enumerate(set(new_clusters_labels))}
#                         new_clusters_labels = [new_cluster_names[cluster] for cluster in new_clusters_labels]
#                         print('Found {0} sub-clusters in cluster {1} with pk of {2}'.format(
#                             len(set(new_clusters_labels)), cluster, pk))
#                         new_clusters = contig_df2cluster_dict(cluster_contig_df, new_clusters_labels)
#                         return [new_clusters, cluster]
#                     pk = int(pk * 0.75)
#                 else:
#                     # print('Failed to find sub-bclusters for {0} with DBSCAN: Could not estimate reachability distance.'.format(shorten_cluster_names(cluster)))
#                     return [{}, cluster]
#         elif cluster_mode == "OPTICS":
#             # print('Trying to find sub-clusters for {0} with DBSCAN.'.format(cluster))
#             # start = timer()
#             # while cluster_contig_df['contig'].size < pk and pk >= min_dims:
#             #     pk = int(pk * 0.75)
#             # if pk < min_dims:
#             #     pk = min_dims
#             # if pk > cluster_contig_df['contig'].size:
#             #     pk = 2
#             # cluster_est = knn_sne_coords(cluster_contig_df, pk)
#             # while not cluster_est and pk >= min_dims:
#             #     pk = int(pk * 0.75)
#             # if pk < min_dims:
#             #     pk = min_dims
#             #     cluster_est = knn_sne_coords(cluster_contig_df, pk)
#             # if cluster_est and cluster_contig_df['contig'].size > pk >= min_dims:
#             #     df_vizbin_coords = cluster_contig_df.loc[:, ['x', 'y']].to_numpy(dtype=np.float64)
#             #     with parallel_backend('threading'):
#             #         dbsc = DBSCAN(eps=cluster_est, min_samples=pk, n_jobs=threads_for_dbscan).fit(df_vizbin_coords)
#             #     new_clusters_labels = dbsc.labels_
#             #     if len(set(new_clusters_labels)) > 1:
#             #         new_cluster_names = {item: cluster + '.' + str(index + 1) for index, item in
#             #                              enumerate(set(new_clusters_labels))}
#             #         new_clusters_labels = [new_cluster_names[cluster] for cluster in new_clusters_labels]
#             #         end = timer()
#             #         print('Found {0} sub-clusters in cluster {1} with DBSCAN and pk of {2} in {3}s.'.format(
#             #             len(set(new_clusters_labels)), cluster, pk, int(end - start)))
#             #         new_clusters = contig_df2cluster_dict(cluster_contig_df, new_clusters_labels)
#             #         return [new_clusters, cluster]
#             # min_samples = int(pow(pk, 2) / 2)
#             #
#             # end = timer()
#             # print('Failed with DBSCAN taking {0}s.Trying to find sub-clusters for {1} using OPTICS with depth.'.format(int(end - start), cluster))
#
#             # print('Trying with OPTICS without depth next.')
#             # while cluster_contig_df['contig'].size < min_samples and min_samples >= min_dims:
#             #     min_samples = int(min_samples * 0.75)
#             # new_clusters_labels = optics_cluster(cluster_contig_df, min_samples=min_samples, include_depth=False,
#             #                                      n_jobs=threads_for_dbscan)[1]
#             # if len(set(new_clusters_labels)) > 1:
#             #     new_cluster_names = {item: cluster + '.' + str(index + 1) for index, item in
#             #                          enumerate(set(new_clusters_labels))}
#             #     new_clusters_labels = [new_cluster_names[cluster] for cluster in new_clusters_labels]
#             #     print('Found {0} sub-clusters in cluster {1} with min samples of {2}'.format(
#             #         len(set(new_clusters_labels)), cluster, min_samples))
#             #     new_clusters = contig_df2cluster_dict(cluster_contig_df, new_clusters_labels)
#             #     return [new_clusters, cluster]
#
#             start = timer()
#             # min_samples = int(pow(pk, 2) / 2)
#             min_samples = pk * 2
#             # Try with automatic pk first
#             min_dims = 4 * 3
#             optics_tries = 0
#             try_interval = 1
#             while len(set(new_clusters_labels)) == 1 and min_dims <= min_samples < cluster_contig_df['contig'].size and optics_tries <= 100:
#                 if optics_tries == try_interval * 10:
#                     try_interval += 1
#                     min_samples = int(min_samples * 0.75)
#                     print('Currently at OPTICS try {0} with min samples of {1} for cluster {2}.'.format(optics_tries,
#                                                                                                                min_samples,
#                                                                                                                shorten_cluster_names(cluster)))
#                 new_clusters_labels = optics_cluster(cluster_contig_df, min_samples=min_samples, include_depth=include_depth,
#                                                      n_jobs=threads_for_dbscan)[1]
#                 if len(set(new_clusters_labels)) > 1:
#                     new_cluster_names = {item: cluster + '.' + str(index + 1) for index, item in
#                                          enumerate(set(new_clusters_labels))}
#                     new_clusters_labels = [new_cluster_names[cluster] for cluster in new_clusters_labels]
#                     end = timer()
#                     # print('Found {0} sub-clusters in cluster {1} with OPTICS and min samples of {2} after {3}s.'.format(
#                     #     len(set(new_clusters_labels)), cluster, min_samples, int(end - start)))
#                     new_clusters = contig_df2cluster_dict(cluster_contig_df, new_clusters_labels)
#                     return [new_clusters, cluster]
#                 optics_tries += 1
#                 min_samples = int(min_samples * 0.75)
#             end = timer()
#             # print('Failed to find sub-bclusters for {0} with OPTICS after {1}s.'.format(shorten_cluster_names(cluster), int(end - start)))
#             return [{}, cluster]
#     else:
#         if clust_dat[5] == 0:
#             # print('Could not calculate purity for cluster {0}. Leaving at 0 and skipping.'.format(shorten_cluster_names(cluster)))
#             return [{}, cluster]
#         elif clust_dat[6] < completeness_threshold:
#             # print('Cluster {0} below completeness threshold with {1}. Skipping.'.format(shorten_cluster_names(cluster),
#             #                                                                             clust_dat[6]))
#             return [{}, cluster]
#         else:
#             print('Cluster {0} meets purity threshold of {1} with {2} and has completeness of {3}.'.format(shorten_cluster_names(cluster),
#                                                                                                            purity_threshold,
#                                                                                                            clust_dat[5],
#                                                                                                            clust_dat[6]))


def uni_dip_depth(clust_dat, cluster, alpha):
    start = timer()
    intervals = [1]
    while len(intervals) == 1 and alpha <= 0.1:
        intervals = UniDip(clust_dat[2], alpha=alpha, ntrials=1000, mrg_dst=1).run()
        if len(intervals) > 1:
            new_clusters = create_new_clusters(cluster, intervals, clust_dat)
            end = timer()
            print('Found {0} depth sub-clusters in cluster {1} with alpha of {2} after {3}s.'.format(
                len(intervals), shorten_cluster_names(cluster), round(alpha, 3), int(end - start)))
            return [new_clusters, cluster]
        alpha += 0.025
    end = timer()
    print('Failed to find depth sub-bclusters for {0} after {1}s.'.format(shorten_cluster_names(cluster),
                                                                          int(end - start)))
    return [{}, cluster]


def dbscan_sub_clusters(cluster_contig_df, cluster, pk, min_dims, threads_for_dbscan, max_tries=10):
    start = timer()

    while cluster_contig_df['contig'].size < pk and pk >= min_dims:
        pk = int(pk * 0.75)
    if pk < min_dims:
        pk = min_dims
    if pk > cluster_contig_df['contig'].size:
        pk = 2

    dbscan_tries = 0
    new_clusters_labels = [1]

    print('Working on {0}.'.format(shorten_cluster_names(cluster)))

    while len(set(new_clusters_labels)) == 1 and min_dims <= pk < cluster_contig_df['contig'].size and dbscan_tries <= max_tries:
        dbscan_tries += 1
        cluster_est = knn_sne_coords(cluster_contig_df, pk)
        while not cluster_est and pk >= min_dims:
            pk = int(pk * 0.75)
        if pk < min_dims:
            pk = min_dims
            cluster_est = knn_sne_coords(cluster_contig_df, pk)
        if cluster_est and cluster_contig_df['contig'].size > pk >= min_dims:
            dims = [dim for dim in ['x', 'y', 'z'] if dim in cluster_contig_df.columns]
            df_sne_coords = cluster_contig_df.loc[:, dims].to_numpy(dtype=np.float64)
            with parallel_backend('threading'):
                dbsc = DBSCAN(eps=cluster_est, min_samples=pk, n_jobs=threads_for_dbscan).fit(df_sne_coords)
            new_clusters_labels = dbsc.labels_
            if len(set(new_clusters_labels)) > 1:
                new_cluster_names = {item: cluster + '.' + str(index + 1) for index, item in
                                     enumerate(set(new_clusters_labels))}
                new_clusters_labels = [new_cluster_names[cluster] for cluster in new_clusters_labels]
                new_clusters = contig_df2cluster_dict(cluster_contig_df, new_clusters_labels)
                end = timer()
                print('Found {0} sub-clusters with DBSCAN in cluster {1} with pk of {2} after {3}s'.format(
                    len(set(new_clusters_labels)), cluster, pk, int(end - start)))
                return [new_clusters, cluster]
            pk = int(pk * 0.75)
        else:
            end = timer()
            print('Failed to find sub-clusters for {0} with DBSCAN after {1}s: Could not estimate reachability distance.'.format(
                shorten_cluster_names(cluster), int(end - start)))
            return [{}, cluster]
    end = timer()
    print(
        'Failed to find sub-clusters for {0} with DBSCAN after {1}s.'.format(
            shorten_cluster_names(cluster), int(end - start)))
    return [{}, cluster]


def optics_sub_clusters(cluster_contig_df, cluster, pk, threads_for_dbscan, include_depth=False, max_tries=50):
    start = timer()
    min_samples = int(pk**2 / 2)
    # min_samples = pk * 2
    # Try with automatic pk first
    dims = [dim for dim in ['x', 'y', 'z'] if dim in cluster_contig_df.columns]
    if include_depth:
        min_dims = 2 * (len(dims) +1)
    else:
        min_dims = 2 * len(dims)
    optics_tries = 0
    try_interval = 1
    new_clusters_labels = [1]

    while len(set(new_clusters_labels)) == 1 and min_dims <= min_samples < cluster_contig_df['contig'].size and optics_tries <= max_tries:
        if optics_tries == try_interval * 10:
            try_interval += 1
            min_samples = int(min_samples * 0.75)
        new_clusters_labels = optics_cluster(cluster_contig_df, min_samples=min_samples, include_depth=include_depth,
                                             n_jobs=threads_for_dbscan)[1]
        if len(set(new_clusters_labels)) > 1:
            new_cluster_names = {item: cluster + '.' + str(index + 1) for index, item in
                                 enumerate(set(new_clusters_labels))}
            new_clusters_labels = [new_cluster_names[cluster] for cluster in new_clusters_labels]
            end = timer()
            print(
                'Found {0} sub-clusters in cluster {1} with OPTICS with include_depth={2} and min samples of {3} after {4}s.'.format(
                    len(set(new_clusters_labels)), cluster, include_depth, min_samples, int(end - start)))
            new_clusters = contig_df2cluster_dict(cluster_contig_df, new_clusters_labels)
            return [new_clusters, cluster]
        optics_tries += 1
        min_samples = int(min_samples * 0.75)
    end = timer()
    print('Failed to find sub-bclusters for {0} with OPTICS with include_depth={1} after {2}s.'.format(shorten_cluster_names(cluster), include_depth, int(end - start)))
    return [{}, cluster]


def hdbscan_sub_clusters(cluster_contig_df, cluster, pk, threads_for_dbscan, depth=False, max_tries=10):
    start = timer()

    dbscan_tries = 0
    new_clusters_labels = [1]

    print('Working on {0}.'.format(shorten_cluster_names(cluster)))

    while len(set(new_clusters_labels)) == 1 and dbscan_tries <= max_tries:
        dbscan_tries += 1

        dims = [dim for dim in ['x', 'y', 'z'] if dim in cluster_contig_df.columns]
        if depth:
            dims.append('depth')
        df_sne_coords = cluster_contig_df.loc[:, dims].to_numpy(dtype=np.float64)
        with parallel_backend('threading'):
            hdbsc = hdbscan.HDBSCAN(core_dist_n_jobs=threads_for_dbscan, min_cluster_size=pk).fit(df_sne_coords)
        new_clusters_labels = hdbsc.labels_
        if len(set(new_clusters_labels)) > 1:
            new_cluster_names = {item: cluster + '.' + str(index + 1) for index, item in
                                 enumerate(set(new_clusters_labels))}
            new_clusters_labels = [new_cluster_names[cluster] for cluster in new_clusters_labels]
            new_clusters = contig_df2cluster_dict(cluster_contig_df, new_clusters_labels)
            end = timer()
            print('Found {0} sub-clusters with HDBSCAN in cluster {1} after {2}s'.format(
                len(set(new_clusters_labels)), cluster, int(end - start)))
            return [new_clusters, cluster]

        else:
            end = timer()
            print('Failed to find sub-clusters for {0} with HDBSCAN after {1}s.'.format(
                shorten_cluster_names(cluster), int(end - start)))
            return [{}, cluster]
    end = timer()
    print(
        'Failed to find sub-clusters for {0} with HDBSCAN after {1}s.'.format(
            shorten_cluster_names(cluster), int(end - start)))
    return [{}, cluster]


def get_sub_clusters2(cluster_dicts, threads_for_dbscan, purity_threshold=0.9, completeness_threshold=0.4, alpha=0.025, pk=None, cluster_mode=None, include_depth=False):
    outputs = []
    for cluster_dict in cluster_dicts:
        cluster = list(cluster_dict.keys())[0]
        # Create dict with just cluster and sort again, to ensure order by depth is
        cluster_dict = sort_cluster_dict_data_by_depth({cluster: cluster_dict[cluster]})
        # All data needed stored in list with following order:
        # 0:contigs, 1:genes, 2:depth, 3:x, 4:y, 5:purity, 6:completeness OR
        # 0:contigs, 1:genes, 2:depth, 3:x, 4:y, 5:z, 6:purity, 7:completeness
        clust_dat = gather_cluster_data(cluster, cluster_dict)
        if len(clust_dat) == 7:
            clust_pur = clust_dat[5]
            clust_comp = clust_dat[6]
        elif len(clust_dat) == 8:
            clust_pur = clust_dat[6]
            clust_comp = clust_dat[7]
        if clust_pur < purity_threshold and isinstance(clust_pur, float) and clust_comp > completeness_threshold:
            print('Cluster {0} below purity of {1} with {2} and above completeness of {3} with {4}. '
                  'Attempting to split.'.format(shorten_cluster_names(cluster), purity_threshold, clust_pur, completeness_threshold, clust_comp))
            uni_dip_out = uni_dip_depth(clust_dat, cluster, alpha)
            if uni_dip_out[0]:
                outputs.append(uni_dip_out)
                continue

            # initialise some stuff
            cluster_contig_df = contig_df_from_cluster_dict(cluster_dict)
            dims = [dim for dim in ['x', 'y', 'z'] if dim in cluster_contig_df.columns]
            min_dims = 2 * len(dims)
            if min_dims > cluster_contig_df['contig'].size:
                min_dims = 2
            if not pk:
                pk = int(np.log(cluster_contig_df['contig'].size))
            if pk < min_dims:
                pk = min_dims

            if cluster_mode == 'DBSCAN' or not cluster_mode:
                print('Trying to find sub-clusters for {0} with DBSCAN.'.format(shorten_cluster_names(cluster)))
                dbscan_out = dbscan_sub_clusters(cluster_contig_df, cluster, pk, min_dims, threads_for_dbscan)
                if dbscan_out[0]:
                    outputs.append(dbscan_out)
                    continue
                outputs.append([{}, cluster])
                continue

            elif cluster_mode == "OPTICS":
                print('Trying to find sub-clusters for {0} with DBSCAN.'.format(cluster))
                dbscan_out = dbscan_sub_clusters(cluster_contig_df, cluster, pk, min_dims, threads_for_dbscan, max_tries=1)
                if dbscan_out[0]:
                    outputs.append(dbscan_out)
                    continue

                print('Trying with OPTICS without depth next for cluster {0}.'.format(cluster))
                optics_out = optics_sub_clusters(cluster_contig_df, cluster, pk, threads_for_dbscan, include_depth=False,
                                                 max_tries=1)
                if optics_out[0]:
                    outputs.append(optics_out)
                    continue

                print('Trying with OPTICS with depth next for cluster {0}.'.format(cluster))
                optics_out = optics_sub_clusters(cluster_contig_df, cluster, pk, threads_for_dbscan, include_depth=True,
                                    max_tries=50)
                if optics_out[0]:
                    outputs.append(optics_out)
                    continue

                outputs.append([{}, cluster])
                continue

            elif cluster_mode == 'HDBSCAN':
                print('Trying to find sub-clusters for {0} with HDBSCAN.'.format(shorten_cluster_names(cluster)))
                hdbscan_out = hdbscan_sub_clusters(cluster_contig_df, cluster, pk, threads_for_dbscan)
                if hdbscan_out[0]:
                    outputs.append(hdbscan_out)
                    continue

                print('Trying with HDBSCAN with depth next for cluster {0}.'.format(cluster))
                hdbscan_out = hdbscan_sub_clusters(cluster_contig_df, cluster, pk, threads_for_dbscan, depth=True)
                if hdbscan_out[0]:
                    outputs.append(hdbscan_out)
                    continue

                outputs.append([{}, cluster])
                continue

        else:
            if clust_pur == 0:
                print('Could not calculate purity for cluster {0}. Leaving at 0 and skipping.'.format(shorten_cluster_names(cluster)))
                outputs.append([{}, cluster])
                continue

            elif clust_comp < completeness_threshold:
                print('Cluster {0} below completeness threshold with {1}. Skipping.'.format(shorten_cluster_names(cluster),
                                                                                            clust_comp))
                outputs.append([{}, cluster])
                continue

            else:
                print('Cluster {0} meets purity threshold of {1} with {2} and has completeness of {3}.'.format(shorten_cluster_names(cluster),
                                                                                                               purity_threshold,
                                                                                                               clust_pur,
                                                                                                               clust_comp))
                outputs.append([cluster_dict, cluster])
                continue
    return outputs


def divide_clusters_by_depth2(ds_clstr_dict, threads, min_purity=0.9, min_completeness=0.4, pk=None, cluster_mode=None, include_depth=False, max_tries=15):
    min_purity = min_purity / 100
    min_completeness = min_completeness / 100
    dict_cp = ds_clstr_dict.copy()
    n_tries = 0

    cluster_list = [[i] + [len(dict_cp[i]['contigs'])] for i in list(dict_cp.keys())]
    cluster_list.sort(key=lambda i: i[1], reverse=True)
    start = timer()
    chunks_to_process = [[] for i in range(threads)]
    for i in cluster_list:
        # chunks_to_process.sort(key=lambda i: len(''.join([contig[1] for contig in i])))
        chunks_to_process.sort(key=lambda i: sum([len(cluster_dict[list(cluster_dict.keys())[0]]['contigs']) for cluster_dict in i]))
        chunks_to_process[0].append({i[0]: dict_cp[i[0]]})
    end = timer()
    print('Created load balanced list in {0}s.'.format(end - start))

    # with Parallel(n_jobs=threads) as parallel:
    inner_max_threads = int(threads/len(chunks_to_process))
    if inner_max_threads < 1:
        inner_max_threads = 1

    previous_dict_keys = list(dict_cp.keys())
    no_progress = False
    split_contigs = []

    with parallel_backend("loky", inner_max_num_threads=inner_max_threads):
        while n_tries <= max_tries and chunks_to_process and not no_progress:
            n_tries += 1
            print('Try: {0}'.format(n_tries))
            sub_clstr_res = Parallel(n_jobs=threads)(delayed(get_sub_clusters2)(cluster_dict_list, inner_max_threads,
                                                                              purity_threshold=min_purity,
                                                                              completeness_threshold=min_completeness,
                                                                              pk=pk, cluster_mode=cluster_mode,
                                                                              include_depth=include_depth) for cluster_dict_list in chunks_to_process)

            # print(sub_clstr_res)

            # for i in dict_cp:
            #     print(i, dict_cp[i].keys())

            for outputs in sub_clstr_res:
                for output in outputs:
                    if output:
                        # print('Removing {0}.'.format(output[1]))
                        if output[1]:
                            dict_cp.pop(output[1])
                            split_contigs.append(output[1])
                        # print('Adding sub-clusters of {0}.'.format(output[1]))
                        dict_cp.update(output[0])
            if sorted(previous_dict_keys) == sorted(list(dict_cp.keys())):
                no_progress = True
                continue
            previous_dict_keys = list(dict_cp.keys())
            # cluster_list = [[i] + [len(dict_cp[i]['contigs'])] for res_list in sub_clstr_res  forif res_list for i in list(e[0].keys())]
            cluster_list = []
            for res_list in sub_clstr_res:
                for res in res_list:
                    if res:
                        for cluster in list(res[0].keys()):
                            cluster_list.append([cluster, len(res[0][cluster]['contigs'])])
            cluster_list.sort(key=lambda i: i[1], reverse=True)
            start = timer()
            # print(cluster_list)
            chunks_to_process = [[] for i in range(threads)]
            for i in cluster_list:
                # chunks_to_process.sort(key=lambda i: len(''.join([contig[1] for contig in i])))
                chunks_to_process.sort(key=lambda i: sum(
                    [len(cluster_dict[list(cluster_dict.keys())[0]]['contigs']) for cluster_dict in i]))
                # print(i[0], dict_cp[i[0]].keys())
                chunks_to_process[0].append({i[0]: dict_cp[i[0]]})
            # for chunk in chunks_to_process:
            #     for some_dict in chunk:
            #         for e in some_dict:
            #             print(e, some_dict[e].keys())
            end = timer()
            print('Created load balanced list in {0}s.'.format(end - start))

    print('Tries until end: {0}'.format(n_tries))
    if dict_cp.keys():
        if 'z' not in list(dict_cp[list(dict_cp.keys())[0]].keys()):
            clust_dat_pur_ind = 5
            clust_dat_comp_ind = 6
        else:
            clust_dat_pur_ind = 6
            clust_dat_comp_ind = 7
    else:
        clust_dat_pur_ind = 6
        clust_dat_comp_ind = 7
    dict_cp = {k: v for k, v in dict_cp.items() if gather_cluster_data(k, dict_cp)[clust_dat_pur_ind] >= min_purity and
               gather_cluster_data(k, dict_cp)[clust_dat_comp_ind] >= min_completeness}
    return dict_cp, split_contigs


# def divide_clusters_by_depth(ds_clstr_dict, threads, min_purity=0.9, min_completeness=0.4, pk=None, cluster_mode=None, include_depth=False):
#     min_purity = min_purity / 100
#     min_completeness = min_completeness / 100
#     dict_cp = ds_clstr_dict.copy()
#     n_tries = 0
#     clusters_to_process = list(dict_cp.keys())
#     # with Parallel(n_jobs=threads) as parallel:
#     inner_max_threads = int(threads/len(clusters_to_process))
#     if inner_max_threads < 1:
#         inner_max_threads = 1
#     with parallel_backend("loky", inner_max_num_threads=inner_max_threads):
#         while n_tries < 15 and clusters_to_process:
#             n_tries += 1
#             print('Try: {0}'.format(n_tries))
#             sub_clstr_res = Parallel(n_jobs=len(clusters_to_process))(delayed(get_sub_clusters)(str(cluster), dict_cp,
#                                                                                                 inner_max_threads,
#                                                                                                 purity_threshold=min_purity,
#                                                                                                 completeness_threshold=min_completeness,
#                                                                                                 pk=pk, cluster_mode=cluster_mode,
#                                                                                                 include_depth=include_depth)
#                                                                       for cluster in clusters_to_process)
#             for output in sub_clstr_res:
#                 if output:
#                     # print('Removing {0}.'.format(output[1]))
#                     dict_cp.pop(output[1])
#                     # print('Adding sub-clusters of {0}.'.format(output[1]))
#                     dict_cp.update(output[0])
#             clusters_to_process = [i for e in sub_clstr_res if e for i in list(e[0].keys())]
#     print('Tries until end: {0}'.format(n_tries))
#     dict_cp = {k: v for k, v in dict_cp.items() if gather_cluster_data(k, dict_cp)[5] >= min_purity and
#                gather_cluster_data(k, dict_cp)[6] >= min_completeness}
#     return dict_cp


def contig_df_from_cluster_dict(cluster_dict):
    clust_contig_df_rows = []
    dims = [dim for dim in ['x', 'y', 'z'] if dim in list(cluster_dict[list(cluster_dict.keys())[0]].keys())]
    # print(dims)
    clust_contig_df_cols = ['contig', 'essential'] + dims + ['depth', 'cluster', 'purity', 'completeness']
    clust_contig_df_ind = []
    clust_contig_df_ind_init = 0
    for i in cluster_dict:
        cluster_essential_genes = [gene for genes in cluster_dict.get(i, {}).get('essential')
                                   for gene in genes.split(',') if not gene == 'non_essential']
        cluster_unique_essential_genes = set(cluster_essential_genes)
        if cluster_essential_genes:
            cluster_purity = round(len(cluster_unique_essential_genes) / len(cluster_essential_genes), 2)
            cluster_completeness = round(len(cluster_unique_essential_genes) / 115, 2)
        else:
            cluster_purity = 0
            cluster_completeness = 0
        for index, contig in enumerate(cluster_dict[i]['contigs']):
            if len(dims) == 2:
                new_row = [contig, cluster_dict[i]['essential'][index], cluster_dict[i]['x'][index],
                           cluster_dict[i]['y'][index], cluster_dict[i]['depth'][index], i, cluster_purity, cluster_completeness]
            elif len(dims) == 3:
                new_row = [contig, cluster_dict[i]['essential'][index], cluster_dict[i]['x'][index],
                           cluster_dict[i]['y'][index], cluster_dict[i]['z'][index], cluster_dict[i]['depth'][index], i, cluster_purity,
                           cluster_completeness]
            clust_contig_df_rows.append(new_row)
            clust_contig_df_ind.append(clust_contig_df_ind_init)
            clust_contig_df_ind_init += 1
    clust_contig_df = pd.DataFrame(clust_contig_df_rows, clust_contig_df_ind,
                                          clust_contig_df_cols)
    return clust_contig_df


def write_scatterplot(df, file_path, hue):
    palette = sns.color_palette("husl", len(set(hue)))
    sorted_set = []
    for i in hue:
        if i not in sorted_set:
            sorted_set.append(i)
    for index, i in enumerate(sorted_set):
        if i == 'N':
            palette[index] = ('silver')
    if 'z' in df.keys():
        df['cluster'] = hue
        palette_dict = {cluster: color for cluster, color in zip(sorted_set, palette)}
        fig = plt.figure()
        ax = Axes3D(fig)
        x = df['x']
        y = df['y']
        z = df['z']
        ax.scatter(x, y, z, s=1.5, c=df['cluster'].map(palette_dict))
        fake_list =[]
        for i, e in zip(sorted_set, palette):
            fake_dot = ax.scatter(df['x'].iloc[0], df['y'].iloc[0], df['z'].iloc[0], s=1.5, color=e)
            fake_list.append(fake_dot)
        ax.legend(fake_list, [str(i) for i in sorted_set], markerscale=2)
        plt.savefig(file_path)
        plt.clf()
    else:
        scatter_plot = sns.scatterplot(data=df, x="x", y="y", hue=hue, sizes=2.5, s=2., palette=palette)
        scatter_plot.legend(fontsize=3, title="Clusters", title_fontsize=4, ncol=1,
                            bbox_to_anchor=(1.01, 1), borderaxespad=0)
        scatter_plot.get_figure().savefig(file_path)
        scatter_plot.get_figure().clf()  # Clear figure


def load_fasta(fasta):
    sequence_dict = {}
    with open(fasta) as f:
        current_header = None
        for line in f:
            # Get header
            if line.startswith('>'):
                # Transform sublists of sequences from previous header to one list.
                if current_header:
                    sequence_dict[current_header] = ''.join(sequence_dict[current_header])
                line = line.strip('\n >')
                sequence_dict[line] = []
                current_header = line
            # Get sequences
            else:
                line = line.strip('\n ')
                sequence_dict[current_header].append(line)
        # Transform sublists of sequences from last header to one list.
        sequence_dict[current_header] = ''.join(sequence_dict[current_header])
    return sequence_dict


def shorten_cluster_names(cluster):
    new_name = []
    sub_clusters = cluster.split('.')
    previous_sub_clstr = None
    for index, sub_clstr in enumerate(sub_clusters):
        if sub_clstr == previous_sub_clstr:
            continue
        previous_sub_clstr = sub_clstr
        counter = 0
        while index+counter+1 <= len(sub_clusters)-1 and sub_clstr == sub_clusters[index+counter+1]:
            counter += 1
        if counter > 0:
            new_name.append(sub_clstr+'x'+str(counter+1))
        else:
            new_name.append(sub_clstr)
    new_name = '.'.join(new_name)
    return new_name


def write_bins(cluster_dict, assembly, min_comp=40, min_pur=90, bin_dir='bins'):
    assembly_dict = load_fasta(assembly)
    # Create bin folder, if it doesnt exist.
    bin_dir = Path(bin_dir)
    bin_dir.mkdir(parents=True, exist_ok=True)
    for cluster in cluster_dict:
        cluster_essential_genes = [gene for genes in cluster_dict.get(cluster, {}).get('essential')
                                   for gene in genes.split(',') if not gene == 'non_essential']
        cluster_unique_ssential_genes = set(cluster_essential_genes)
        if cluster_essential_genes:
            cluster_purity = round(len(cluster_unique_ssential_genes) / len(cluster_essential_genes) * 100, 0)
            cluster_completeness = round(len(cluster_unique_ssential_genes) / 115 * 100, 0)
        else:
            cluster_purity = 0
            cluster_completeness = 0
        if cluster_purity >= min_pur and cluster_completeness >= min_comp:
            new_cluster_name = shorten_cluster_names(cluster)
            bin_name = '_'.join(['binny'] + ['C'+str(cluster_completeness)] + ['P'+str(cluster_purity)] + [new_cluster_name])
            bin_file_name = bin_name + '.fasta'
            bin_out_path = bin_dir / bin_file_name
            with open(bin_out_path, 'w') as out_file:
                for contig in cluster_dict[cluster]['contigs']:
                    out_file.write('>' + contig + '\n' + assembly_dict.get(contig) + '\n')

########################################################################################################################

def rec_comp(kmer):
    base_pairs = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    rc_kmer = "".join(base_pairs.get(base, base) for base in reversed(kmer))
    return rc_kmer


def kmer_is_canonical(kmer):
    rc_kmer = rec_comp(kmer)
    if kmer <= rc_kmer:
        canonical_kmer = kmer
        is_canonical = True
    else:
        canonical_kmer = rc_kmer
        is_canonical = False
    return canonical_kmer, is_canonical
    # return is_canonical


def kmer_counter(kmer, sequence):
    d = {i: ['A', 'T', 'G', 'C'] for i in range(kmer)}
    kmer_list = [''.join(combination) for combination in itertools.product(*[d[k] for k in sorted(d.keys())])]
    kmer_list_can = [i for i in kmer_list if kmer_is_canonical(i)[1]]
    kmer_count = [[kl_km, len([i for i in re.finditer(r'(?=({0}))'.format(kl_km), sequence)]) +
                   len([i for i in re.finditer(r'(?=({0}))'.format(rec_comp(kl_km)), sequence)])] for kl_km in kmer_list_can]
    return kmer_count


def frequencyCount(string, substr):
   count = 0
   pos = 0
   while(True):
       pos = string.find(substr , pos)
       if pos > -1:
           count = count + 1
           pos += 1
       else:
           break
   return count


def kmer_counter2(kmer_list_can, sequence):
    # kmer_count = [[kl_km, len([i for i in re.finditer(r'(?=({0}))'.format(kl_km), sequence)]) +
    #                len([i for i in re.finditer(r'(?=({0}))'.format(rec_comp(kl_km)), sequence)])] if kl_km != rec_comp(kl_km)
    #               else [kl_km, len([i for i in re.finditer(r'(?=({0}))'.format(kl_km), sequence)])] for kl_km in kmer_list_can]
    # kmer_count = [[kl_km, frequencyCount(sequence, kl_km) + frequencyCount(sequence, rec_comp(kl_km))]
    #               if kl_km != rec_comp(kl_km) else [kl_km, frequencyCount(sequence, kl_km)] for kl_km in kmer_list_can]
    kmer_count = [(frequencyCount(sequence, kl_km) + frequencyCount(sequence, rec_comp(kl_km))) / len(sequence)
                  if kl_km != rec_comp(kl_km) else frequencyCount(sequence, kl_km) / len(sequence) for kl_km in kmer_list_can]
    return kmer_count


def kmer_counter2_bloom(ksize, kmer_list_can, sequence):
    target_table_size = 5e8
    num_tables = 4
    bloomfilter = khmer.Nodetable(ksize, target_table_size, num_tables)
    bloomfilter.consume(sequence)
    # kmer_count = [bloomfilter.get(kl_km) / len(sequence)
    #               if bloomfilter.get(kl_km) > 0 and bloomfilter.get(rec_comp(kl_km)) == 0
    #               else bloomfilter.get(rec_comp(kl_km)) / len(sequence)
    #               if bloomfilter.get(kl_km) == 0 and bloomfilter.get(rec_comp(kl_km)) > 0
    #               else 0 for kl_km in kmer_list_can]
    # kmer_count = [[kl_km, bloomfilter.get(kl_km)] for kl_km in kmer_list_can]
    kmer_count = [bloomfilter.get(kl_km) / len(sequence) for kl_km in kmer_list_can]
    return kmer_count


def get_contig_kmer_freq2(ksize, can_k_mers, assembly_chunk):
    chunk_list = []
    for contig in assembly_chunk:
        start = timer()
        if ksize < 7:
            kmer_count_list = kmer_counter2(can_k_mers, contig[1])
        else:
            kmer_count_list = kmer_counter2_bloom(ksize, can_k_mers, contig[1])
        end = timer()
        if end - start >= 10:
            print('Counting k-mers in {0} took longer than 10s with {1}s.'.format(contig[0], end - start))
        # print('Counted k-mers in {0} in {1}s.'.format(contig[0], end - start))
        contig_kfreq = [contig[0]] + kmer_count_list
        # contig_kfreq = kmer_count_list
        # contig_kfreq_med = np.median(contig_kfreq[1:])
        # contig_kfreq_97quant = np.quantile(contig_kfreq[1:], q=0.97)
        # contig_kfreq = [contig[0]] + [count_graph.get(i) / kmers_found if count_graph.get(i) < contig_kfreq_97quant
        #                               else contig_kfreq_med / kmers_found for i in range(nkmers)]
        chunk_list.append(contig_kfreq)
    return chunk_list


def get_contig_kmer_matrix2(contig_list, ksize_list, n_jobs=1):
    contig_kmer_freq_matrix = []

    contig_list = [i + [len(i[1])] for i in contig_list]
    contig_list.sort(key=lambda i: i[2], reverse=True)
    start = timer()
    chunks_to_process = [[] for i in range(n_jobs)]
    for i in contig_list:
        # chunks_to_process.sort(key=lambda i: len(''.join([contig[1] for contig in i])))
        chunks_to_process.sort(key=lambda i: sum([contig[2] for contig in i]))
        chunks_to_process[0].append(i)
    end = timer()
    print('Created load balanced list in {0}s.'.format(end - start))
    # Try to free mem
    del contig_list
    n_rounds = 0
    print('k-mer sizes to count: {0}'.format(', '.join([str(i) for i in ksize_list])))
    for ksize in ksize_list:
        d, kmer_list, kmer_list_can = {}, [], []

        d = {i: ['A', 'T', 'G', 'C'] for i in range(ksize)}
        kmer_list = [''.join(combination) for combination in itertools.product(*[d[k] for k in sorted(d.keys())])]
        kmer_list_can = sorted(list({i for i in kmer_list if kmer_is_canonical(i)[1]}))

        # contig_list.sort(key=lambda i: len(i[1]), reverse=True)
        # for i, e in zip(enumerate(contig_list),
        #                 [e for i in range(-(-len(contig_list) // n_jobs)) for e in list(range(n_jobs))]):
        #     chunks_to_process[e].append(contig_list[i[0]])
        if n_rounds == 0:
            contig_kmer_freq_matrix.append(['contig']+kmer_list_can)
        else:
            contig_kmer_freq_matrix[0] = contig_kmer_freq_matrix[0] + kmer_list_can
        with parallel_backend("loky"):
            contig_kmer_freq_matrix_chunks = Parallel(n_jobs=n_jobs)(delayed(get_contig_kmer_freq2)(ksize, kmer_list_can, chunks) for chunks in chunks_to_process)
        contig_counter = 0
        for chunk in contig_kmer_freq_matrix_chunks:
            for contig_freq in chunk:
                contig_counter += 1
                if n_rounds == 0:
                    contig_kmer_freq_matrix.append(contig_freq)
                else:
                    if contig_kmer_freq_matrix[contig_counter][0] != contig_freq[0]:
                        print(contig_kmer_freq_matrix[contig_counter][0],contig_freq[0])
                        raise Exception
                    contig_kmer_freq_matrix[contig_counter] = contig_kmer_freq_matrix[contig_counter] + contig_freq[1:]
        print('Finished counting k-mer frequencies for size {0}'.format(ksize))
        n_rounds += 1
    return contig_kmer_freq_matrix


def split(a, n):
    k, m = divmod(len(a), n)
    return (a[i * k + min(i, m):(i + 1) * k + min(i + 1, m)] for i in range(n))


def get_perp(n_data_points):
    return int(np.log(n_data_points**(np.log(n_data_points)/np.log(2))))
    # return int(np.log(n_data_points**np.log(n_data_points)))


def gff2low_comp_feature_dict(annotation_file):
    contig_dict = {}
    # Get all info for each line
    with open(annotation_file, 'r') as af:
        for line in af:
            line = line.strip(' \t\n').split('\t')
            contig, feature_start, feature_end, attributes = line[0], line[3], line[4], line[8].split(';')
            # Check if feature has no start and/or stop position or length = 0
            if not feature_start or not feature_end or int(feature_start)-int(feature_end) == 0:
                continue
            # Only add line if feature was annotated as essential gene
            for attribute in attributes:
                attribute_name, attribute_id = attribute.split('=')[0], attribute.split('=')[1]
                if attribute_name == 'product' and 'ribosomal RNA' in attribute_id:
                    if not contig_dict.get(contig):
                        contig_dict[contig] = [(int(feature_start), int(feature_end))]
                    else:
                        contig_dict[contig].append((int(feature_start), int(feature_end)))
                elif attribute_name == 'rpt_family':
                    if not contig_dict.get(contig):
                        contig_dict[contig] = [(int(feature_start), int(feature_end))]
                    else:
                        contig_dict[contig].append((int(feature_start), int(feature_end)))
    # # Create list of lists to join dict items as string and build data frame
    # annot_cont_line_list = [[contig, ','.join(contig_dict.get(contig))] for contig in contig_dict]
    # annot_cont_df = pd.DataFrame(annot_cont_line_list, columns=['contig', 'essential'])
    return contig_dict


def mask_rep_featrues(contig_rep_feature_dict, contig_list):
    for index, contig in enumerate(contig_list):
        if contig_rep_feature_dict.get(contig[0]):
            regions = contig_rep_feature_dict.get(contig[0])
            # print(contig[0], len(contig[1]), regions)
            for region in regions:
                if region[0] - 1 <= 0:
                    start_region = 0
                    start2feature = 0
                else:
                    start_region = region[0] - 1
                    start2feature = region[0] - 2
                contig_list[index][1] = contig[1][:start2feature]\
                                    + contig[1][start_region:(region[1] - 1)].lower()\
                                    + contig[1][region[1]:]
                # contig[1] = contig_rep_masked
                # print(index, contig[0], region, contig[1][(region[0] - 5):(region[0]+5)], '...', contig[1][(region[1] - 5):(region[1]+5)])
                # print(contig)
    # return contig_list


def binny_iterate(contig_data_df, threads, min_purity, min_completeness):
    leftovers_df = contig_data_df.copy()
    n_iterations = 0
    n_new_clusters = 1
    good_clusters = {}
    while n_iterations <= 10 and n_new_clusters > 0:
        n_iterations += 1
        print('Running iteration: {0}.'.format(n_iterations))
        # Initial clustering
        # start = timer()
        init_clust_dict, labels = run_initial_scan(leftovers_df, 'HDBSCAN', threads, include_depth=False)
        # Write first scan to file
        inititial_clust_df = cluster_df_from_dict(init_clust_dict)
        inititial_clust_df.to_csv('intermediary/iteration_{0}_initial_contigs2clusters.tsv'.format(n_iterations), sep='\t', index=False)
        # end = timer()
        # print(end - start)

        # Plot initial final clustering
        write_scatterplot(leftovers_df, 'intermediary/iteration_{0}_initial_scatter_plot.pdf'.format(n_iterations),
                          labels)

        # Find sub-clusters
        print('Attempting to divide leftover clusters by depth with UniDip and/or sub-clusters by subsequent HDBSCAN runs.')
        new_clust_dict, split_clusters = divide_clusters_by_depth2(init_clust_dict, threads, int(min_purity),
                                                                   int(min_completeness), cluster_mode='HDBSCAN',
                                                                   include_depth=True, max_tries=10)
        new_clust_dict = {'R{0}.'.format(n_iterations) + k: v for k, v in new_clust_dict.items()}
        n_new_clusters = len(set(new_clust_dict.keys()))
        # round_leftovers_dict = {k: v for k, v in init_clust_dict.items() if k not in split_clusters}
        # iteration_clust_dict = {**new_clust_dict, **round_leftovers_dict}
        iteration_clust_dict = {**new_clust_dict, **init_clust_dict}
        # iteration_clust_dict = {**new_clust_dict}

        good_clusters.update(new_clust_dict)
        good_clust_contig_df = contig_df_from_cluster_dict(good_clusters)

        iteration_clust_df = cluster_df_from_dict(iteration_clust_dict)
        iteration_clust_df.to_csv('intermediary/iteration_{0}_contigs2clusters.tsv'.format(n_iterations), sep='\t', index=False)

        # iteration_clust_contig_dft = new_clust_contig_df.merge(contig_data_df, how='outer', copy=False)

        iteration_clust_contig_df = contig_df_from_cluster_dict(iteration_clust_dict)

        # Plot sub-clustering
        conditions = [(iteration_clust_contig_df['purity'] >= min_purity / 100) & (
                iteration_clust_contig_df['completeness'] >= min_completeness / 100),
                      (iteration_clust_contig_df['purity'] < min_purity / 100) | (
                              iteration_clust_contig_df['completeness'] < min_completeness / 100)]
        values = [iteration_clust_contig_df['cluster'], 'N']
        iteration_clust_contig_df['above_thresh'] = np.select(conditions, values)
        iteration_clust_contig_df = contig_data_df.merge(iteration_clust_contig_df, how='outer', on='contig',
                                                         suffixes=(None, '_y'))
        iteration_clust_contig_df['above_thresh'] = iteration_clust_contig_df['above_thresh'].fillna('N')

        iteration_clust_contig_df = iteration_clust_contig_df[
            ~((iteration_clust_contig_df['contig'].duplicated(keep=False))
              & (iteration_clust_contig_df['above_thresh'] == 'N'))]

        write_scatterplot(iteration_clust_contig_df,
                          'intermediary/iteration_{0}_final_scatter_plot.pdf'.format(n_iterations),
                          iteration_clust_contig_df['above_thresh'])
        leftovers_df = iteration_clust_contig_df[iteration_clust_contig_df['above_thresh'] == 'N']
        leftovers_df = leftovers_df[~leftovers_df['contig'].isin(good_clust_contig_df['contig'].tolist())]
    return good_clusters, init_clust_dict
