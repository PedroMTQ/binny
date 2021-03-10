import pandas as pd
import seaborn as sns
from pathlib import Path

proj_dir = Path('/Users/oskar.hickl/binny_bench/binner_data/checkm/checkm_results_001')
binny = list(proj_dir.glob('binny_out/WA-?_v001/checkm/storage/bin_stats_ext.tsv'))
concoct = list(proj_dir.glob('concoct_out/WA-?_v001/checkm/storage/bin_stats_ext.tsv'))
maxbin = list(proj_dir.glob('maxbin_out/WA-?_v001/checkm/storage/bin_stats_ext.tsv'))
metabat = list(proj_dir.glob('metabat_out/WA-?_v001/checkm/storage/bin_stats_ext.tsv'))
files_list = [binny, concoct, maxbin, metabat]

for files in files_list:
    binner = str(files[0]).split('/')[7].split('_')[0]
    binner_good_bins = 0
    binner_conts = []
    binner_completenesses = []
    for file in files:
        with open(file, 'r') as f:
            for line in f:
                line = line.strip('\n \t')
                line_bin = line.split('\t')[0]
                line_data = line.split('\t')[1]
                line_data = line_data.strip('{}').split(',')
                for data in line_data:
                    if 'Completeness' in data:
                        line_comp = float(data.split(': ')[1])
                        binner_completenesses.append(line_comp)
                    elif 'Contamination' in data:
                        line_cont = float(data.split(': ')[1])
                        binner_conts.append(line_cont)
                if line_comp >= 40 and line_cont <= 10:
                    binner_good_bins += 1
                    # print(line_bin, line_comp, line_cont)
    avg_cont = round(sum(binner_conts) / len(binner_conts), 1)
    avg_comp = round(sum(binner_completenesses) / len(binner_completenesses), 1)
    avg_good_bins = int(binner_good_bins / 4)
    print(binner, avg_good_bins, avg_cont, avg_comp)




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
        print('Running initial scan with DBSCAN and min samples of: {0}'.format(str(pk)))
        # Run parallelized dbscan
        first_clust_dict, labels = hdbscan_cluster(contig_data_df, pk=pk, n_jobs=dbscan_threads,
                                                  include_depth=include_depth)

    return first_clust_dict, labels


def hdbscan_cluster(contig_data_df, pk=None, include_depth=False, n_jobs=1):
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
        hdbsc = hdbscan.HDBSCAN(min_cluster_size=15, core_dist_n_jobs=n_jobs).fit(dim_df)
    cluster_labels = hdbsc.labels_
    cluster_dict = contig_df2cluster_dict(contig_data_df, cluster_labels)
    return cluster_dict, cluster_labels









