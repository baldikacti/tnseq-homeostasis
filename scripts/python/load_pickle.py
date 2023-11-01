import glob
import pickle
import numpy as np
import pandas as pd
import multiprocessing

unique_counts_caulo = ['data/peter/peter_UC_processed.csv']


for i in range(len(unique_counts_caulo)):
	data_file = unique_counts_caulo[i]
	df_data = pd.read_csv(data_file)
	data_SM = df_data.iloc[:,7:]
	data_SM = np.asarray(data_SM)
	data = data_SM
	print(data.shape)

	nthread = 4
	n_pickle = nthread
	# n_pickle = 6

	ls_ess = []
	ls_non_ess = []
	for i in range(n_pickle):
		with open(data_file + 'essential_' + str(i) + '.pkl', 'rb') as f:
			temp_ess = pickle.load(f)
			for j in range(len(temp_ess)):
				if temp_ess[j] == 'failed' or temp_ess[j] < 0:
					temp_ess[j] = 0


		with open(data_file + 'non_essential_' + str(i) + '.pkl', 'rb') as f:
			temp_non_ess = pickle.load(f)
			for j in range(len(temp_non_ess)):
				if temp_non_ess[j] == 'failed':
					temp_non_ess[j] = 500

		ls_ess.append(temp_ess)
		ls_non_ess.append(temp_non_ess)

	ls_ess_flat = [item for sublist in ls_ess for item in sublist]
	thresold_ess = np.asarray(ls_ess_flat).reshape(-1, 1)

	ls_non_ess_flat = [item for sublist in ls_non_ess for item in sublist]
	thresold_non_ess = np.asarray(ls_non_ess_flat).reshape(-1, 1)

	thresold = np.column_stack((thresold_ess, thresold_non_ess))
	thresold = thresold.astype('int')
	print(thresold)
	print(thresold.shape)
    
	with open(data_file + '_' + 'threshold.npy', 'wb') as f:
		np.save(f, thresold)


	assert ('failed' in thresold) == False
	assert np.shape(thresold) == (len(data[0]), 2)
