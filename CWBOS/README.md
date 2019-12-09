# CW-BOS imaging of high intensity ultrasonic pressure fields.


In this study, there are three parts to quantitatively map focuced-ultrasound (FUS) pressure fields based on CW-BOS imaging.

1. [CW-BOS simulations](https://github.com/wgrissom/zebrography/blob/master/CWBOS/CWBOS_simulations): Perform numerical FUS beam simulations and generate the training data. 
  *demo_Simulations.m  
    Workflow to generate simulated data, perform SVD on the dictionary, and construct training set.
  *wave_prop_simu.m
    Generate simulated pressure of different levels and transducers. 
  *forward_model_dxdz.m
    Generate displacements in both two dismentions on the background pattern.
  *forward_model_proj.m
    Generate the projected pressure waveform.
  *accum_d.m
    Caculate displacements of each slice, need to be called in forward_model_dxdz.m

2. [CW-BOS acquistion and hardware](https://github.com/wgrissom/zebrography/tree/master/CWBOS/CWBOS_acquistion): Acquire FUS-photo by CW-BOS system. 

3. [Reconstruction](https://github.com/wgrissom/zebrography/tree/master/CWBOS/recon): Train deep neural network, process the acquired acutal photo by CW-BOS system and reconstruct root-mean-square(RMS) projected pressure by pre-trained neural network. 
