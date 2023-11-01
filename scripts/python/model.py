import numpy as np
import numpy.random as rn
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import time
from tqdm import tqdm_notebook
import glob
import scipy.stats as ss

from path import Path

import pyximport
pyximport.install()

from apf.models.gap_hdpseq import GaP_HDP_Seq

datafile_ls = glob.glob("data/peter/*.csv")
# threshold_ls = glob.glob("threshold/*.npy")

for datafile in datafile_ls:
	print(datafile)

	### process data ###
	df_data = pd.read_csv(datafile)
	df_data = df_data.sort_values(by=['strain', 'condition']).reset_index(drop=True)
	df_data.to_csv('sort_data.csv', index=False)
	# df_data[df_data>30] = 30
	# df_data.to_csv('y_data.csv', index=False)
	print(df_data)
	data = df_data.iloc[:,7:]
	# data[data>50] = 50
	data = np.asarray(data)
	print(data)
	datafile = datafile[11:][:-4] + '_'
	print(datafile)
	# break
	# temp_path = datafile[16:][:-14] + datafile[16:][-4:]
	temp = np.load("data/peter/peter_UC_processed.csv_threshold.npy")
	# temp = np.load('data/threshold_npy/thresold.npy')
	# temp = np.repeat(np.array([[1,50]]), data.shape[1], axis=0)
	# temp = np.array([[1,50], [1, 50], [1, 50], [1, 10], [1, 10],
	# 				 [1, 30], [1,100], [1,100], [1,30], [1,10]])
	print(temp.shape)
	for i in range(len(temp)):
		if temp[i][0] == 0:
			temp[i][0] = 1
		# temp[i] = [0,10]

		if abs(temp[i][0] - temp[i][1]) <= 5 or temp[i][0] - temp[i][1] > 0:
	    	# temp[i][0] = 1
			temp[i][1] = max(temp[i][0], temp[i][1])
		if temp[i][1] - temp[i][0] < 10:
			temp[i][1] = temp[i][0] + 10
		if np.mean(data[:,i]) < 2:
			temp[i] = [1, 10]

	np.save('threshold.npy', temp)
	assert np.load("threshold.npy").shape == (data.shape[1], 2)

	
	### Model ###
	data_SM = data
	n_subpops = 1000
	L = data.shape[1]
	# n_biopsies_I = np.asarray([2,4])
	n_biopsies_I = np.asarray(df_data.groupby(['strain','condition']).size()) ### this sort by alphabeta
	assert sum(n_biopsies_I) == data.shape[0]

	seed = np.random.randint(0,9999)
	infer_model = GaP_HDP_Seq(data_SM = data_SM,
	                          n_loci=L,
	                          n_biopsies_I=list(n_biopsies_I),
	                          n_subpops=n_subpops,
	                          shp_gam=1.,
	                          rte_gam=1.,
	                          shp_beta=1.,
	                          rte_beta=1.,
	                          seed=seed,
	                          n_threads=60)

	# burnin
	path_name = datafile + '_psamples'
	out_dir = Path(path_name)
	out_dir.makedirs_p()

	n_itns = 10
	n_total_itns = 10000
	n_epochs = int(n_total_itns / n_itns)

	# fix_state = {'gam': 1e-3,
	#              'beta': 1e3}
	start = time.time()
	for epoch in range(n_epochs+1):
	    infer_model.fit(data_SM,
	                    n_itns=0 if epoch==0 else n_itns,
	                    initialize=(epoch==0),
	                    verbose=n_itns)
	#                     fix_state=fix_state)

	    state = dict(infer_model.get_state())
	    np.savez_compressed(out_dir.joinpath('%d_state.npz' % infer_model.total_itns), **state)

	end = time.time()
	print('running time is %d second' % (end-start))
	print('seed is %i' % seed)

	
	### Results ###
	samples = glob.glob(datafile + '_psamples/*state*')
	print('%d posterior samples' % len(samples))

	H = []
	G = []
	Gi = []
	Gij = []

	total_sample = sum(n_biopsies_I)
	# threshold = 0.0005
	n_psample = len(samples)

	for s in range(len(samples)):
	    sample_file = np.load(samples[s])  # load in sample
	    # this is the H matrix (using paper notation), sorry for different notation
	    h = sample_file['Z_LK'].T
	#    print(Z_LK.shape)          # this is a loci x component matrix

	    Theta_K = sample_file['Theta_K']   # unnormalized G distribution
	    # unnormalized per-individual G distributions
	    Theta_IK = sample_file['Theta_IK']
	    # unnormalized per-sample G distributions
	    Theta_SK = sample_file['Theta_SK']

	    # normalized G distribution
	    g = Theta_K / Theta_K.sum()
	    # normalized per-individual G distributions
	    gi = Theta_IK / np.sum(Theta_IK, axis=1, keepdims=True)
	    # normalized per-sample G distributions
	    gij = Theta_SK / np.sum(Theta_SK, axis=1, keepdims=True)

	    H.append(h)
	    G.append(g)
	    Gi.append(gi)
	    Gij.append(gij)

	flattern = [str(s.tolist()) for item in H for s in item]

	glist = [s for item in G for s in item]
	# glist = [s for item in G[-1] for s in item]



	gilist = []
	for i in range(len(n_biopsies_I)):
	    flattern_gi = [s for item in Gi for s in item[i]]
	    gilist.append(flattern_gi)


	gijlist = []
	for i in range(total_sample):
	    flattern_gij = [s for item in Gij for s in item[i]]
	    gijlist.append(flattern_gij)


	df = pd.DataFrame({'H': flattern})

	df['g0'] = glist

	for i in range(len(gilist)):
	    df['gi_'+str(i+1)] = gilist[i]

	for i in range(len(gijlist)):
	    df['gij_'+str(i+1)] = gijlist[i]

	df_temp = df.groupby(['H']).agg('sum').reset_index()

	df_temp.to_csv('temp.csv', index=False)

	# threshold = 0.00099
	# df_p = df_temp[(df_temp[df_temp.columns[-(total_sample):]]
	#                 > n_psample*threshold).any(axis=1)]

	# df_p[df_p.columns[1:len(df_p.columns)]] = df_p[df_p.columns[1:len(
	#     df_p.columns)]].apply(lambda x: x / sum(x))

	# df_p = df_p.sort_values('g0', ascending=False)

	# df_p.to_csv('inference_threshold.csv', index=False)
	path_name_2 = 'results/' + datafile
	out_dir = Path(path_name_2)
	out_dir.makedirs_p()

	df_temp[df_temp.columns[1:len(df_temp.columns)]] = df_temp[df_temp.columns[1:len(
    df_temp.columns)]].apply(lambda x: x / sum(x))

	# df_p = df_temp.nlargest(5, 'g0')
	threshold = 0.0005
	# df_p = df_temp[(df_temp[df_temp.columns[-(total_sample):]]
	# > threshold).any(axis=1)]
	df_p = df_temp[df_temp['g0']> threshold]
	df_p[df_p.columns[1:len(df_p.columns)]] = df_p[df_p.columns[1:len(
    df_p.columns)]].apply(lambda x: x / sum(x))

	df_p = df_p.sort_values('g0', ascending=False)

	row = ['H'] + ['Top_Level'] + \
		  list(df_data.groupby(['strain', 'condition']).groups.keys()) + \
		  list(df_data.groupby(['strain', 'condition', 'slevel', 'rep']).groups.keys())
	
	row_condition = pd.DataFrame([row], columns=list(df_p.columns))

	df_annote =  pd.concat([row_condition, df_p])

	df_annote.to_csv('results/%s/%s_inference_oak_1.csv' % (datafile, datafile), index=False)
	
	locus_att = pd.read_csv('locus_attribute.tab', sep='\t')
	locus_att = locus_att[['locus_tag', 'Name', 'description']].drop_duplicates(subset=['locus_tag'])
	locus_tag = df_data.columns[7:]
	df_locus_tag = pd.DataFrame(locus_tag)
	df_locus_tag.columns = ['locus_tag']
	df_att = df_locus_tag.merge(locus_att, on='locus_tag', how='left')
	
	n_cluster = len(df_p)
	H_matrix = np.zeros((n_cluster,L), dtype='int')
	H_matrix_str = np.asarray(df_p['H'])
	for i in range(n_cluster):
		H_matrix_str[i] = H_matrix_str[i].replace(", ", "").replace("[", "").replace("]", "")

	for i in range(n_cluster):
		for j in range(L):
			if H_matrix_str[i][j] == "0":
				H_matrix[i,j] = int(0)
			if H_matrix_str[i][j] == "1":
				H_matrix[i,j] = int(1)

	df_H = pd.DataFrame(H_matrix.T)
	df_H["locus_tag"] = df_att['locus_tag']
	df_H['Gene_Name'] = df_att['Name']
	df_H['Description'] = df_att['description']

	cols = list(df_H.columns)
	cols = [cols[-3]] + [cols[-2]] + [cols[-1]] + cols[:-3]
	df_H = df_H[cols]


	df_H.to_csv('results/%s/%s_H_matrix_oak_1.csv' %(datafile, datafile), index=False)
	# df_H.to_csv('results/' + datafile + '/' + datafile + '_H_matrix.csv', index=False)



