#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Training data writer.

This scirpt is used to compress *.mat data into *.hdf5.

Created on Mon Dec  2 12:44:51 2019

@author: huiwenluo
"""
import h5py
import numpy as np
Basepath = './'
dict1 = h5py.File(Basepath+'/dict0910svd116.mat','r')


#training set
xtrain = np.transpose(np.float32(dict1['hisdic']))
V_red = np.float32(dict1['V_red'])
with h5py.File(Basepath+'/dict_train0910svd116.hdf5', 'w') as f:
    dset = f.create_dataset("train", data=xtrain,compression='gzip',compression_opts=9)

#Target set: RMS projected pressure.
ytrain = np.transpose(np.float32(dict1['rmsproj']),[1,0])
with h5py.File(Basepath+'/dict_target0910svd116.hdf5', 'w') as f:
    dset = f.create_dataset("rmsproj", data=ytrain,compression='gzip',compression_opts=9)
    
#Put projection matrix in info.
with h5py.File(Basepath+'/dict_info0910svd116.hdf5', 'w') as f:
    dset = f.create_dataset("V_red", data=V_red.T,compression='gzip',compression_opts=9)
