    #!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Train a simple deep neural network with multiple hidden layers.

Input: projected histograms based on the compressed SVD space
Output: RMS projected pressure


"""
#%%
from __future__ import print_function
from keras.layers import Dense
from keras.optimizers import RMSprop
from keras.callbacks import ModelCheckpoint, LearningRateScheduler
from keras.callbacks import ReduceLROnPlateau
from keras.models import Sequential
from keras.layers.convolutional import regularizers
import numpy as np
import os
import h5py
from random import shuffle
import svd_generator

dict1 =h5py.File('/scratch/luoh3/dict_train0910svd116.hdf5','r')
x_train= np.float32(dict1['train'])
dict2 =h5py.File('dict_target0910svd116.hdf5','r')
rmsprojtrain= np.float32(dict2['rmsproj'])
print(np.size(rmsprojtrain))
dict3 =h5py.File('dict_info0910svd116.hdf5','r')
V_red= np.float32(dict3['V_red'])
y_train = rmsprojtrain

print(y_train.max())
print(y_train.min())
print(rmsprojtrain.max())
print(rmsprojtrain.min())

x_train = np.reshape(x_train,[x_train.shape[0],54,54,1])
x_train = np.transpose(x_train,(0,2,1,3))
# Input image dimensions.
input_shape = x_train.shape[1:]
print('x_train shape:', x_train.shape)
print(x_train.shape[0], 'train samples')
print('y_train shape:', y_train.shape)
y_train = (y_train-y_train.mean())/y_train.std()
ind_list = [i for i in range(np.size(x_train,axis = 0))]
shuffle(ind_list)
x_train = x_train[ind_list,:]
y_train = y_train[ind_list,:]

def lr_schedule(epoch):
    """Learning Rate Schedule

    Learning rate is scheduled to be reduced after 80, 120, 160, 180 epochs.
    Called automatically every epoch as part of callbacks during training.

    # Arguments
        epoch (int): The number of epochs

    # Returns
        lr (float32): learning rate
    """
    lr = 5e-5
    print('Learning rate: ', lr)
    return lr
    
name = 'model' 
n1 = 2e-5
n2 = 2e-4
dim = V_red.shape[1]
model = Sequential()
model.add(Dense(dim,kernel_regularizer=regularizers.l1_l2(l1 = n1, l2 = n2),kernel_initializer='he_normal',input_shape=(dim,),activation = 'tanh'))
model.add(Dense(dim,kernel_regularizer=regularizers.l1_l2(l1 =  0, l2 = n2),kernel_initializer='he_normal',activation = 'tanh'))
model.add(Dense(dim,kernel_regularizer=regularizers.l1_l2(l1 =  0, l2 = n2),kernel_initializer='he_normal',activation = 'tanh'))
model.add(Dense(dim,kernel_regularizer=regularizers.l1_l2(l1 =  0, l2 = n2),kernel_initializer='he_normal',activation = 'tanh'))
model.add(Dense(1,kernel_initializer='he_normal',kernel_regularizer=regularizers.l1_l2(l1 =  0, l2 = n2),activation = 'linear'))


model.summary() 
# print the model summary

# determine Loss function and Optimizer
model.compile(loss='mean_squared_error',
              optimizer=RMSprop(lr=lr_schedule(0)),
              metrics=['mean_squared_error'])

# Prepare model model saving directory.
save_dir = os.path.join(os.getcwd(), 'model_dnn')
model_name = 'zebra_%s_model.{epoch:03d}.h5'
if not os.path.isdir(save_dir):
    os.makedirs(save_dir)
filepath = os.path.join(save_dir, model_name)
 
# Prepare callbacks for model saving and for learning rate adjustment.
checkpoint = ModelCheckpoint(filepath=filepath,
                             monitor='loss',
                             verbose=1,
                             save_best_only=True)

lr_scheduler = LearningRateScheduler(lr_schedule)

lr_reducer = ReduceLROnPlateau(factor=np.sqrt(0.1),
                               cooldown=0,
                               patience=5,
                               min_lr=0.5e-6)

callbacks = [checkpoint, lr_reducer, lr_scheduler]

batch_size =1024
epochs = 20

#This is a modifed data generator to project vectorized histograms into the compressed SVD subspace.

batch_size =1024
num_classes = 1
epochs = 20
import svd_generator
datagen = svd_generator.ImageDataGenerator(
    rotation_range=30,
    validation_split = 0,
    height_shift_range = 0,
    width_shift_range = 0,
  #  brightness_range=[0.7,1.3]
    V_red = V_red
    )
datagen.fit(x_train)
split = int(len(x_train))
train_generator = datagen.flow(x_train[:split,:,:,:], y_train[:split],batch_size = batch_size,shuffle = False)
history = model.fit_generator(generator=train_generator,
                 verbose=2,epochs = epochs, steps_per_epoch = 400,
                    callbacks=callbacks)
