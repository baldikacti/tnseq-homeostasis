#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul  5 15:06:04 2021

@author: shaihe
"""

import glob
import pickle
import pymc3 as pm
import pandas as pd
import numpy as np
import theano
import theano.tensor as tt

import multiprocessing
from multiprocessing import Process


def sampling(data_file, block_i, block_len):
    thresold_ess = []
    thresold_non_ess = []
    df_data = pd.read_csv(data_file)
    data = df_data.iloc[:,5:]
    df_data = df_data.iloc[:, 0:11]
    data = np.asarray(data)
    print(data.shape)
    for i in range(block_i[0],block_i[-1]+1):
        X = data[:,i]
        try:
            mixture_fix = pm.Model()
            with mixture_fix as model:
                p = pm.Dirichlet("p", a=np.array([1.0, 1.0]), shape=2)
                non_essential_mean = pm.Normal("non_essential_mean", mu=100, sigma=10)
                essential_mean = pm.Normal("essential_mean", mu=0, sigma=1)
                means = tt.as_tensor([essential_mean, non_essential_mean])
                sds = tt.as_tensor([1, 10])
            #    sds = tt.stack([1,1])
            #    means = tt.stack([0, non_essential_mean])
                cat = pm.Categorical("category", p=p, shape=len(X))
                # sd = pm.Uniform("sd", lower=0, upper=20)
                obs = pm.Normal('obs', mu=means[cat], sigma=sds[cat], observed=X)
                
            with model:
                tr = pm.sample(5000, cores=1, chains=1, progressbar = False)
            ess_score = np.mean(tr['essential_mean'])
            non_ess_score = np.mean(tr['non_essential_mean'])
            thresold_ess.append(ess_score)
            thresold_non_ess.append(non_ess_score)

        except Exception:
            thresold_ess.append('failed')
            thresold_non_ess.append('failed')

    print("essential is", thresold_ess)
    print("non essential is", thresold_non_ess)

    with open(data_file +'essential_'+ str(block_i[0]//block_len) +'.pkl', 'wb') as f:
       pickle.dump(thresold_ess, f)

    with open(data_file + 'non_essential_' + str(block_i[0]//block_len) + '.pkl', 'wb') as f:
       pickle.dump(thresold_non_ess, f)

if __name__ == '__main__':

    unique_counts_caulo = ['data/adam/Caulo_unique_processed.csv']
    print(len(unique_counts_caulo))

    for i in range(len(unique_counts_caulo)):
        data_file = unique_counts_caulo[i]
        df_data = pd.read_csv(data_file)
        df_data = df_data.iloc[:, 0:11]
        data = df_data.iloc[:,5:]
        data = np.asarray(data)
        
        print(data_file)
        print(data)
        nrow = len(data)
        ncol = len(data[0])
    
        nthread = 6
        loci_ls = list(range(ncol))
        a, b = divmod(ncol, nthread-1)
        ls = []
        for i in range(0, len(loci_ls), a):
            if i == len(loci_ls) - b:
                ls.append(loci_ls[i:])
                break

            else:
                ls.append(loci_ls[i:i+a])
#         print(ls)
#         print(len(ls))
        assert len(ls) == min(len(ls), nthread)

#         print(ls)
        block_len = len(ls[0])
        
        proc_ls = []
        for i in range(min(nthread, len(ls))):
            proc = Process(target=sampling, args=([data_file, ls[i], block_len]))
            proc_ls.append(proc)


        for item in proc_ls:
            item.start()

        for item in proc_ls:
            item.join()
       