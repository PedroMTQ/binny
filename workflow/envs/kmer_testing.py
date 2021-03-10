contig_list.sort(key=lambda i: i[0])
kfreq_array.sort(key=lambda i: i[0])
# kfreq_array_multi.sort(key=lambda i: i[0])
ksize = 4
nkmers = 4**ksize
tablesize = nkmers + 10
cg = khmer.Countgraph(ksize, tablesize, 1)
cg.set_use_bigcount(True)

myind = 0
myc = contig_list[myind][1]
cg.consume(myc)

myk = kfreq_array[myind]

print(kfreq_array[-1])
print(kfreq_array[myind])
print(myk)
myl = [[cg.reverse_hash(i), cg.get(i)] for i in range(nkmers)]
myl.sort(key=lambda i: i[0])

khmer_out = [contig_list[myind][0]] + myl
print(khmer_out)

print(khmer_out == myc)

for i in myk[1:]:
  if i not in khmer_out and not [rec_comp(i[0]), i[1]] in khmer_out:
    print(i)

############
#Some other stuff
# super_contig = []
# for i in contig_list[1:]:
#     if sum([len(i) for i in super_contig]) < 10000000:
#         super_contig.append(i[1])
# super_contig = ''.join(super_contig)
# len(super_contig)
# # f_chunk = [['super_contig', super_contig]]
#
# fcl = get_contig_kmer_freq_test(6, f_chunk)
# print(sorted(fcl[0][1:], reverse=True)[1:20])
#
# bad_c = 0
# # for cfreq, contig in zip(kfreq_array, contig_list[1:]):
# for cfreq in kfreq_array:
#     # if 65535 in cfreq[1:]:
#     #     bad_c += 1
#         # print('Contig {0}, length {1} is all fucked up.'.format(cfreq[0], len(contig[1])))
#     # print(sorted(set(cfreq[1:]), reverse=True)[1000:1100])
#     print(len(cfreq[1:]))
# print(bad_c)
###########
# Old/unused function etc
# start = timer()
# kmer_counter2(9, mys)
# end = timer()
# print(end - start)

# def kmer_counter_alt(kmer, sequence):
#     d = {i: ['A', 'T', 'G', 'C'] for i in range(kmer)}
#     kmer_dict = {''.join(combination): 0 for combination in itertools.product(*[d[k] for k in sorted(d.keys())])}
#     kmer_list_can = sorted([i for i in kmer_dict.keys() if kmer_is_canonical(i)[1]])
#     for kmer, count in kmer_dict.items():
#         kmer_dict[kmer_is_canonical(kmer)[0]] += len([i for i in re.finditer(r'(?=({0}))'.format(kmer), sequence)])
#     kmer_dict = {k: v for k,v in kmer_dict.items() if k in kmer_list_can}
#     return kmer_dict

# def get_contig_kmer_freq(ksize, assembly_chunk):
#     chunk_list = []
#     nkmers = 4 ** ksize
#     tablesize = nkmers + 10
#     if ksize > 7:
#         target_table_size = 5e8
#         num_tables = 4
#         count_graph = khmer.Nodetable(ksize, target_table_size, num_tables)
#         for contig in assembly_chunk:
#             count_graph.consume(contig[1])
#             kmers_found = sum([count_graph.get(i) for i in range(nkmers)])
#             contig_kfreq = [contig[0]] + [count_graph.get(i) / kmers_found for i in range(nkmers)]
#             chunk_list.append(contig_kfreq)
#     else:
#         count_graph = khmer.Countgraph(ksize, tablesize, 1)
#         count_graph.set_use_bigcount(True)
#         for contig in assembly_chunk:
#             count_graph.consume(contig[1])
#             kmers_found = sum([count_graph.get(i) for i in range(nkmers)])
#             contig_kfreq = [contig[0]] + [count_graph.get(i) / kmers_found for i in range(nkmers)]
#             chunk_list.append(contig_kfreq)
#     return chunk_list


# def get_contig_kmer_freq_test(ksize, assembly_chunk):
#     chunk_list = []
#     nkmers = 4 ** ksize
#     tablesize = nkmers + 10
#     # if ksize > 7:
#     #     target_table_size = 5e8
#     #     num_tables = 4
#     #     count_graph = khmer.Nodetable(ksize, target_table_size, num_tables)
#     #     for contig in assembly_chunk:
#     #         count_graph.consume(contig[1])
#     #         kmers_found = sum([count_graph.get(i) for i in range(nkmers)])
#     #         contig_kfreq = [contig[0]] + [count_graph.get(i) / kmers_found for i in range(nkmers)]
#     #         chunk_list.append(contig_kfreq)
#     # else:
#     for contig in assembly_chunk:
#         count_graph = khmer.Countgraph(ksize, tablesize, 1)
#         count_graph.set_use_bigcount(True)
#         count_graph.consume(contig[1])
#         kmers_found = sum([count_graph.get(i) for i in range(nkmers)])
#         contig_kfreq = [contig[0]] + [count_graph.get(i) for i in range(nkmers)]
#         # contig_kfreq_med = np.median(contig_kfreq[1:])
#         # contig_kfreq_97quant = np.quantile(contig_kfreq[1:], q=0.97)
#         # contig_kfreq = [contig[0]] + [count_graph.get(i) / kmers_found if count_graph.get(i) < contig_kfreq_97quant
#         #                               else contig_kfreq_med / kmers_found for i in range(nkmers)]
#         chunk_list.append(contig_kfreq)
#     return chunk_list


# def get_contig_kmer_matrix(contig_list, ksize, n_jobs=1):
#     contig_kmer_freq_matrix = []
#     chunks_to_process = [[] for i in range(n_jobs)]
#     contig_list = [i + [len(i[1])] for i in contig_list]
#     contig_list.sort(key=lambda i: i[2], reverse=True)
#     start = timer()
#     for i in contig_list:
#         # chunks_to_process.sort(key=lambda i: len(''.join([contig[1] for contig in i])))
#         chunks_to_process.sort(key=lambda i: sum([contig[2] for contig in i]))
#         chunks_to_process[0].append(i)
#     end = timer()
#     print('Created load balanced list in {0}s.'.format(end - start))
#     # Try to free mem
#     del contig_list
#     # for i in contig_list:
#     #     chunks_to_process.sort(key=lambda i: len(''.join([contig[1] for contig in i])))
#     #     chunks_to_process[0].append(i)
#     nkmers = 4 ** ksize
#     tablesize = nkmers + 10
#     count_graph = khmer.Countgraph(ksize, tablesize, 1)
#     contig_kmer_freq_matrix.append(['contig']+[count_graph.reverse_hash(i) for i in range(nkmers)])
#     with parallel_backend("loky"):
#         contig_kmer_freq_matrix_chunks = Parallel(n_jobs=n_jobs)(delayed(get_contig_kmer_freq_test)(ksize, chunks) for chunks in chunks_to_process)
#     for chunk in contig_kmer_freq_matrix_chunks:
#         for contig_freq in chunk:
#             contig_kmer_freq_matrix.append(contig_freq)
#     return contig_kmer_freq_matrix

#######

# start = timer()
# kfreq_array_1 = get_contig_kmer_matrix2(contig_list, [2], threads)
# end = timer()
# print(end - start)
#
# start = timer()
# kfreq_array_2 = get_contig_kmer_matrix2(contig_list, [4], threads)
# end = timer()
# print(end - start)
#
# start = timer()
# kfreq_array_3 = get_contig_kmer_matrix2(contig_list, 5, threads)
# end = timer()
# print(end - start)

# kfreq_array = [kfreq_array_1[index] + c_k_freq[1:] for index, c_k_freq in enumerate(kfreq_array_2)]
# kfreq_array = [kfreq_array[index] + c_k_freq[1:] for index, c_k_freq in enumerate(kfreq_array_3)]

###

# kfreq_array = kfreq_array_1


# kfreq_array = kfreq_array_multi


# #####
# start = timer()
# kfreq_array2 = get_contig_kmer_matrix(contig_list, 2, threads)
# end = timer()
# print(end - start)
# #####

# print(len(kfreq_array))
# print(kfreq_array[1][:20])
# print(kfreq_array[-1][:20])

# print(X_std_scale[0][:10])
# print(X_std_scale[-1][:10])

# transformer = RobustScaler().fit(X)
# X_scaled = transformer.transform(X)
# print(X_rob_scale[0][:10])
# print(X_rob_scale[-1][:10])

####

# if sum(pca.explained_variance_ratio_[:5]) > 0.9:
#     X_pca = X_scaled
# else:
#     X_pca = transformer.transform(X_scaled)

# X_pca = transformer.transform(X_scaled)

# fig = plt.figure()
# df = pd.DataFrame(X_pca, columns=['x', 'y'])
#     df['contig'] = X_contigs
#     sns.scatterplot(data=df, x="x", y="y", sizes=0.4, s=4)
# plt.show()
