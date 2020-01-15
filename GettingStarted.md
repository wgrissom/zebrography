
Getting Started
==
These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. See deployment for notes on how to deploy the project on a live system.

#### All the example data are published in [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3601557.svg)](https://doi.org/10.5281/zenodo.3601557).

### Prerequisites
----

#### Software
* MATLAB
* Python (Tensorflow, Keras and [Keras processing](https://github.com/keras-team/keras-preprocessing))
* [Pythonista](http://omz-software.com/pythonista/) installed on the iPad.

#### Hardware Setup
* An ultra-clear rimless water tank
* An acoustic absorber (Aptflex F48, Precision Acoustic Ltd, UK)
* A mega-pixel digital singe-lens reflex camera
* An iPad or other tablets.
* Focused-ultrasound (FUS) transducer
* Waveform generator
* Amplifier 
* An Arduino board, a wired DSLR camera shutter and an analog switch

### CW-BOS Photo Acquisition Process
---
1. Fill the water tank with degassed deionized water.
2. Connect an Arduino Board to the external port of waveform generator and the analog switch that is connected to the wired switch (see the shutter design in "/CWBOS_acquisition/modified_camera_shutter_design.zip").
3. Connect the wave generator to the amplifier, and connect the amplifier to the focused-ultrasound transducer. Connect the DSLR camera to the experimental computer via USB and open the software to control the camera remotely (e.g.: EOS Utility Software for Canon DSLR camera). 
4. Compile the code "/CWBOS_acquisition/ShutterController_Arduino.ino" on the Arduino board. 
5. Run "BOSTomoDisplay_app_iPad.py" on the Pythonista installed on the iPad,  and keep the IP address shown on the APP in mind!
6. Change the IP address in the "/CWBOS_acquisition/testacq.m' and set other parameters in the script.
5. Run "BOSTomoDisplay_app_iPad.py" on the Pythonista installed on the iPad,  and keep the IP address shown on the APP in mind!
6. Change the IP address in the "/CWBOS_acquistion/testacq.m' and set other parameters in the script.
7. Set the parameters of waveform generator and the camera parameters.
8. Put the black cloth over the water tank and the camera setup to suppress the ambient light.
10. Open the amplifier and then open the output of the waveform generator
11. Run  "/CWBOS_acquistion/testacq.m” on the experiment computer and in the meantime photos can be saved on the experiment computer. 



### Software Execution
---
#### Training set generation
1. Simulate spatially- and temporally-resolved steady-state FUS pressure fields with nonlinearity using "/CWBOS_simulations/wave_prop_simu.m". An example of simulated data is stored in "p0_152500.mat".

```Matlab
p0 = 152500; % Transmitter pressure amplitude  
f0 = 1.1e6; % [Hz] Center frequency of focused-ultrasound  
a = 31.6e-3; % [m] Source radius.  
f_num = 2; % f-number of focused-ultrasound.  
[apaz_sv,dX,dY,dZ,nX,nY,z_sv] = wave_prop_simu(p0,a,f0,f_num);  
```
###### Tips: Save simulated pressure data of each amplitude for convenience in future. 

2.  Calculate projected pressure and physical displacements in meters by "/CW_simulations/forward_model.m". 
```Matlab
Zd = 17.061e-2/2; % Distance between iPad screen and the middle of FUS beam.
[dxreal,dzreal,proj] = forward_model(apaz_sv,dX,dY,dZ,nX,nY,z_sv,Zd);
```
###### "dxreal" and "dzreal" are displacements in x- and z-dimension, respectively. "proj" is projected pressure waveform. 

3. Calculate the point spread function from the actual photo taken by camera. 
You can decide the size of histograms on your own. We decide it and pin width by counting pixels from the actual photo. One example is below:

```Matlab
nblock = 54; % Size of blocks, even number to be consistent with iPad and camera
npin = 6;% pin: 1 pixel, spacing: 8 pixels
blocks = ones(nblock,nblock);
blocks((nblock-npin+2)/2:(nblock+npin)/2,(nblock-npin+2)/2:(nblock+npin)/2) = 0;
block1 = blocks;
load('seg.mat')
block2 = seg.unblur(:,:,8,13:14);
block2 = block2(:,:,:); % collapse third and fourth dimensions
ft1 = fftshift(fftshift(fft2(block1),1),2);
ft2 = fftshift(fftshift(fft2(block2),1),2);
F2 = ft2(:);
F1 = repmat([spdiag(ft1(:)) ones(nblock*nblock,1)],[size(block2,3) 1]);
lambda = 1;
x = (F1'*F1 + lambda*eye(nblock*nblock+1))\(F1'*F2);
win = fspecial('gaussian',54,17);
win = win ./ max(win(:));
xfun = x(1:end-1).*win(:);
figure; imagesc(abs(reshape(xfun,[54,54])))
X = real(ift2(reshape(x(1:end-1),[nblock nblock])));
fits = reshape(F1(:,1:end-1)*xfun(:),[nblock nblock size(block2,3)]);
fits1 = (ifft2((fits(:,:,1))));
```

4. Generate simulated histograms. Simulated histograms of p0 = 152500 and f_num = 2 are given in "his_152500_2.mat".
```Matlab
ds = 25.4/264*1e-3*1668/834/8; %% width of each pixel
P0 = 152500;
nT = 40;
n = 4;
nP0 = length(P0);
nblock = 54; % Size of blocks, even number to be consistent with iPad and camera
npin = 6;% pin: 1 pixel, spacing: 8 pixels
blocks = ones(nblock,nblock);
blocks((nblock-npin+2)/2:(nblock+npin)/2,(nblock-npin+2)/2:(nblock+npin)/2) = 0;
f_num = 2;
for pp = 1:nP0
load(['displacement_',num2str(P0(pp)),'_',num2str(f_num),'.mat']);
locx = round(dxreal/ds);
locz = round(dzreal/ds);
nX = size(locx,2);
nZ = size(locz,3);
his = zeros(nblock,nblock,nX,nZ);
for ii = 1+n/2:nX-n/2-1
    for jj = 1+n/2:nZ-n/2-1
    histmp = zeros(nblock,nblock);
        for kk = 1:nT
        x = round(locx(kk,ii,jj)); z = round(locz(kk,ii,jj));
        histmp = histmp +  circshift(blocks,[x,z]);
        end
    histmp = histmp/nT;
    %-------------------------------------------------------------%
    %Apply the point spread funtion for the simple simulated FUS
    %histograms. xfun is the frequency transform of the calculated
    %point spread function by last section.
    ft1 = fftshift(fftshift(fft2(histmp),1),2);
    F1 = [spdiag(ft1(:)) ones(nblock*nblock,1)];
    fits = reshape(F1*xfun(:),[nblock nblock]);
    fits1 = abs(ifft2(fits(:,:,1)));
    %--------------------------------------------------------------%
    his(:,:,ii,jj) = fits1;
    end
end
his = single(his(:,:,3:end-3,3:end-3));
nx = size(his,3); nz = size(his,4);
%Vectorize rectangular histograms. 
his = his(:,:,:); 
his = permute(his(:,:,:),[3 1 2]);
his = his(:,:);
save(['his_',num2str(P0(pp)),'_',num2str(f_num),'.mat'],'his','nx','nz');
end
```

5. Concatenate histograms, perform the singular value decomposition (SVD) and save projection matrix "V_red". A simple example only including data of p0 = 152500 and f_num = 2 is stored in "dict.mat".

```Matlab
P0 = 152500; %% we saved each simulated dataset by p0 (transmitter pressure) in simulations. 
nP0 = length(P0);
hisdic = [];
rmsproj = [];
fnum = [1,2];
for ff = 1:length(fnum)
    for pp = 1:nP0
        load(['proj_',num2str(P0(pp)),'_',num2str(f_num),'.mat']);
        load(['his_',num2str(P0(pp)),'_',num2str(f_num),'.mat']);
        proj = squeeze(proj(:,3:end-3,3:end-3));
        nx = size(proj,2);
        nz = size(proj,3);
        proj = proj(:,:);
        rms = squeeze(sqrt(mean(proj.^2,1)))';
        histest = his-his(1,:);
        ind = find(sqrt(mean(histest.^2,2))~=0);
        %each image is normalized by subtracting its own mean and being divided
        %by its own standard deviation
        %"ind" is used to delete the duplicated entries.
        hisdic = cat(1,hisdic,(his(ind,:)-mean(his(ind,:),2))./std(his(ind,:),[],2));
        rmsproj = cat(1,rmsproj,rms(ind));
    end
end
hisdic = cat(1,hisdic,(his(1,:)-mean(his(1,:),2))./std(his(1,:),[],2));
rmsproj = cat(1,rmsproj,0);
[~,S,V] = svd(hisdic,'econ');
s_v = diag(S);
e = cumsum(s_v.^2)./sum(s_v.^2);
nDictSpace = sum(e<=1-1e-5);
V_red = V(:,1:nDictSpace);
save(['dict.mat'],'hisdic','rmsproj','V','S','V_red','-v7.3')
```

###### a) "/CW_simulations/demo_simulations.m" includes the whole process to acquire the final training set;
###### b) Python script "/recon/demo_traniningdata_writer.py" can be used to covert and large "*mat" file to "*.hdf5".
7. Run the Python script "/recon/svd_trainDNN.py" to train the neural network.
8. Run the Matlab script "/recon/process_photo.m" to process the actual photos that you acquire in the experiments and save the set of histograms.
9. Run the python script "/recon/demo_predict.py" to reconstruct the root-mean-square projected pressure maps.


