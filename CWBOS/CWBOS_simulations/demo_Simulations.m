%| Demonstration of wave simulations, SVD and making training sets using
%| a small dataset.
%| Copyright 2019, Huiwen Luo and William Grissom, Vanderbilt University.

addpath('./')
addpath('../util')
p0 = 152500; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulate FUS pressure fields using a modified angular specturm approch
% with non-linearity.
% It takes several hours by current spatial and temporal resolution, and
% FOV. User can set resolutions and FOV on their own. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a=31.6e-3; % source radius, we are simulation H101 transducer here. 
f0=1.16e6; % center frequency of pulse
f_num = 2;
is_save = false;
[apaz_sv,dX,dY,dZ,nX,nY,z_sv] = wave_prop_simu(p0,a,f0,f_num,is_save);

%% 
% Calculate projected pressure and displacements of pixels in the background
% pattern using simulated pressure datasets above. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Zd = 17.061e-2/2; % The distance between iPad screen and the middle of FUS beam.
[apaznew,dxreal,dzreal] = forward_model_dxdz(apaz_sv,dX,dY,dZ,nX,nY,z_sv,Zd);
proj = forward_model_proj(apaz_sv,dY,dZ,nX,nY,z_sv);
save(['./data_eg/dataparams_',num2str(P0(pp)),'_',num2str(f_num),'.mat'],'dX','dY','dZ','nX','nZ','nY','p0','-v7.3');
%save(['pressure_',num2str(p0),'_',num2str(f_num),'.mat'],'apaznew','-v7.3');
save(['./data_eg/displacement_',num2str(p0),'_',num2str(f_num),'.mat'],'dxreal','dzreal','-v7.3');
%% 
% To make simulated histograms closer to the reality, the point spread
% funtion needs to be calculated from the actual photo by the least-square
% fitings. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nblock = 54; % Size of blocks, even number to be consistent with ipad and camera
npin = 6;% pin: 1 pixel, spacing: 8 pixels
blocks = ones(nblock,nblock);
PSF = fspecial('average',3);
blocks((nblock-npin+2)/2:(nblock+npin)/2,(nblock-npin+2)/2:(nblock+npin)/2) = 0;
block1 = blocks;
load('./data_eg/seg.mat')%200mv_1_5.mat')
block2 = seg.unblur(:,:,8,13:14);
block2 = block2(:,:,:); % collapse third and fourth dimensions
ft1 = fftshift(fftshift(fft2(block1),1),2);
ft2 = fftshift(fftshift(fft2(block2),1),2);
F2 = ft2(:);
F1 = repmat([spdiag(ft1(:)) ones(nblock*nblock,1)],[size(block2,3) 1]);
lambda = 1;
x = (F1'*F1 + lambda*eye(nblock*nblock+1))\(F1'*F2);
xfun = x; %%FFT of point spread function 

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Now generate training data of simulated histograms!
% Apply the calculated displacements above into non-FUS histograms, and
% convolve with the point spread function.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ds = 25.4/264*1e-3*1668/834/6; %% width of each pixel
P0 = 152500;
nT = 40;
n = 4;
nP0 = length(P0);
nblock = 54; % Size of blocks, even number to be consistent with ipad and camera
npin = 6;% pin: 1 pixel, spacing: 8 pixels
blocks = ones(nblock,nblock);
blocks((nblock-npin+2)/2:(nblock+npin)/2,(nblock-npin+2)/2:(nblock+npin)/2) = 0;
f_num = 2;
for pp = 1:nP0
    load(['./data_eg/displacement_',num2str(P0(pp)),'_',num2str(f_num),'.mat']);
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

    save(['./data_eg/his_',num2str(P0(pp)),'_',num2str(f_num),'.mat'],'his');
end
%% Make a dictionary comprising of vectorized histograms and peform SVD on the dictionary.

%%Expand intial histograms to SVD space 
%%Load simualted histograms and projected pressure and make the dictionary.
P0 = 152500; %% we saved each simulated dataset by p0(transmitter pressure) in simualtions. 

%For 1.16MHz H101 transducer, P0 = [325:150:3200]*100, f_num = [1,2]
%For 2.25MHz transducer, P0 = [325:150:2800]*100, f_num = [1,2,3];
nP0 = length(P0);
hisdic = [];
for pp = 1:nP0
%     load(['./data_eg/proj_',num2str(P0(pp)),'_',num2str(f_num),'.mat']);
%     nx = size(proj,2);
%     nz = size(proj,3);
    load(['./data_eg/his_',num2str(P0(pp)),'_',num2str(f_num),'.mat']);

    histest = his-his(1,:);
    ind = find(sqrt(mean(histest.^2,2))~=0);
    %each image is normalized by subtracting its own mean and being divided
    %by its own standard deviation
    %"ind" is used to delete the duplicated entries.
    hisdic = cat(1,hisdic,(his(ind,:)-mean(his(ind,:),2))./std(his(ind,:),[],2));
end

%Add a non-FUS histogram entries into the dictionary.
hisdic = cat(1,hisdic,(his(1,:)-mean(his(1,:),2))./std(his(1,:),[],2));

%%Perform SVD, truncate SVD space into a compressed subspace.
[~,S,V] = svd(hisdic,'econ');
s_v = diag(S);
e = cumsum(s_v.^2)./sum(s_v.^2);
nDictSpace = sum(e<=1-1e-5);
V_red = V(:,1:nDictSpace); %V_red here is the truncated right singular matrix. 


%Then we project histograms into the truncated SVD space. 
P0 = 152500;
nP0 = length(P0);
rmsproj = [];
train = [];
for pp = 1:nP0
    %     data = load(['dict_',num2str(P0(pp)),'.mat']);
    load(['./date_demo/proj_',num2str(P0(pp)),'_',num2str(f_num),'.mat']);
    load(['./date_demo/his_',num2str(P0(pp)),'_',num2str(f_num),'.mat']);
    proj = squeeze(proj(:,3:end-3,3:end-3));
    nx = size(proj,2);
    nz = size(proj,3);
    proj = proj(:,:);
    rms = squeeze(sqrt(mean(proj.^2,1)))';
    his = (his(:,:)-mean(his(:,:),2))./std(his(:,:),[],2);
    train = cat(1,train,his*V_red);
    rmsproj = cat(1,rmsproj,rms(:));
end

save(['./date_demo/dict.mat'],'train','rmsproj','S','V','V_red','-v7.3');




