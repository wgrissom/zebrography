# CW-BOS imaging of high intensity ultrasonic pressure fields.


In this study, there are three parts to quantitatively map focuced-ultrasound (FUS) pressure fields based on CW-BOS imaging.

1. [CW-BOS simulations](https://github.com/wgrissom/zebrography/blob/master/CWBOS/CWBOS_simulations): Perform numerical FUS beam simulations and generate the training data. 
  - demo_Simulations.m  
    * Workflow to generate simulated data, perform SVD on the dictionary, and construct training set.
  - wave_prop_simu.m 
    * Generate simulated pressure of different levels and transducers. 
    "ablvec.p", "march_asr.p"  "precalculate_abl", "precalculate_ad.p" and "precalculate_mas.p" need to be called.
  - forward_model_dxdz.m
    * Generate displacements in both two dismentions on the background pattern.
  - forward_model_proj.m
    * Generate the projected pressure waveform.
  - accum_d.m
    * Caculate displacements of each slice, need to be called in forward_model_dxdz.m

2. [CW-BOS acquistion and hardware](https://github.com/wgrissom/zebrography/tree/master/CWBOS/CWBOS_acquistion): Acquire FUS-photo by CW-BOS system. 
  - testacq
    * Workflow to acquire FUS and non-FUS photos.
    * Need to call "BOSTomoController.py" to communicate with iPad. 
  - BOSTomoController.py 
    * To be ran on the experiment computer to switch the background pattern.
  - BOSTomoDisplay_app_IPad.py 
    * To be ran on an app named [Pythonista](http://omz-software.com/pythonista/) of iPad to display background patterns and IP address of iPad.
  - ShutterController_Arduino.ino
    * Code to be uploaded to the Arduino board, which allows Arduino to control waveform generator and DLSR camera.
  - Modified_camera_shutter_design.zip
    * Modified camera shutter with an analog switch. 
     

3. [Reconstruction](https://github.com/wgrissom/zebrography/tree/master/CWBOS/recon): Train deep neural network, process the acquired acutal photo by CW-BOS system and reconstruct root-mean-square(RMS) projected pressure by pre-trained neural network.
  - demo_trainingdata_writer.py
    * Write training data with HDF5 format.
  - demo_svd_trainDNN.py
    * Train a multi-layer deep neural network.
  - process_photo.m
    * Process actual photos acquired by DSLR camera. 
  - demo_predict.py
    * Reconstruct RMS projected pressure from actual photos. 
  - "model116.h5" and "model225.h5" are pre-trained model for two transducers (1.16MHz and 2.25MHz). 
