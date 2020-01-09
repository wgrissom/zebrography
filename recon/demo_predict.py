
"""
Reconstruct RMS projected pressure map by bringing segmented histograms from the actual photo into
pre-trained neural network.

Created in 2019

@author: huiwenluo, Vanderbilt University, 2019
"""
#%%   
import scipy.io as sio 
import numpy as np
import h5py
import os
from keras.models import load_model  
import svd_generator

# Load the training set to fit the date generator. 
dict1 =h5py.File('/dict_targe1213svd116.hdf5','r')
ytrain= np.float32(dict1['rmsproj'])
dict3 =h5py.File('/dict_info1213svd116.hdf5','r')
V_red= np.float32(dict3['V_red'])
dict1 =h5py.File('/dict_train1213svd116.hdf5','r')
xtrain= np.float32(dict1['train'])
xtrain=np.reshape(xtrain,[xtrain.shape[0],54,54,1])
xtrain = np.transpose(xtrain,(0,2,1,3))
testdatagen = svd_generator.ImageDataGenerator(
    V_red = V_red
   )
testdatagen.fit(xtrain)

# Load the pre-trained model and predict RMS projected pressure. 
# Before bringing data into the neural network, remember that histograms need to be normalized by subtracting the mean and being divided by the STD individually.
loadmodel = load_model('model116.h5') # Load the pre-trained neural network. 
repnum = 1 # Number of repetitions
for pn in range(repnum):
    dict1 = sio.loadmat('tilehis200mv'+str(pn+1)+'.mat') # load the histograms segmented from the actual photos. 
    nx = int(dict1['nX'])
    nz = int(dict1['nZ'])
    histblur = np.float32(dict1['tilehisblur'])
    histblur = np.transpose(histblur,(2,3,0,1))
    histblur = np.reshape(histblur,(nx*nz,2916)).T
    histblur = (histblur-np.mean(histblur,axis = 0))/np.std(histblur,axis = 0) # normalize each histograms by subtracting its mean and dividing it by its standard deviation.
    histblur = np.reshape(histblur.T,(nx*nz,54,54,1))
    test_generator = testdatagen.flow(histblur,shuffle = False,batch_size = histblur.shape[0]) # project the histograms to the SVD subspace. 
    rmsblur = loadmodel.predict_generator(test_generator,steps = 1,verbose = 0)*ytrain.std()+ytrain.mean() # predict the corresponding root-mean-square projected pressure.
    rmsblur =  np.reshape(rmsblur,[nx,nz]) # reshape the vector to reconstructed root-mean-square projected map. 
    sio.savemat('predict'+str(pn+1)+'.mat',{'rmsblur':rmsblur,'nx':nx,'nz':nz,'histblur':histblur}) # save the map as *mat format file.

  
