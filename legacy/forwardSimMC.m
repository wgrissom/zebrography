%% 
% parameters
f = 1.16e6; % Hz, US frequency
omega = 2*pi*f; % rad/sec, US frequency
c = 1540; % m/s, US speed in water
lambda = c/f; % m, US wavelength in water
k = 2*pi/lambda; % wavenumber
n0 = 1.3325; % base index of refraction
Zd = 0.09; % m, distance from center of US focus to background
dn0dp = (1.3475-1.3325)/((1100-1)*100000); % dindex/dPascal, from Waxler paper, read from Figure 2 at 24.8 deg C
pmax = 1.65E5; % Pascal; peak pressure, at 200mV
%Texp = 1/30; % exposure time
%dt = 1/f/19.9999; % 20 samples per US period
Nt = 100; % time samples
%t = 0:dt:Texp - dt; % time vector
dt = 1/f/Nt;
t = 0:dt:1/f - dt;
pix = 0.178/1920; %pixel size = length / no. of pixel

% load a field profile; define spatial vector
load FOCUSout.mat %%%load H101_Sim_1p1MHz; 
% this field is rotationally symmetric but we need to shift it slightly
p = real(P3d);
p = p./max(abs(p(:)))*pmax; %%%pmap_single = pmap_single./max(pmap_single(:))*pmax;
%%%p = squeeze(pmap_single(:,:,:)); % remove 3 points to center focus for tomo recon
[Nx,Ny,Nt] = size(p);
% interpolate to 10x finer res in y
% for ii = 1:Ny
%     for jj = 1:Nz
%         pmap_highres(:,ii,jj)
%     end
% end
dx = pix; % resolution in m, 0.064 m FOV comes from Charles' email
dy = pix;
dz = dy;
%x = 0:dx:(Nx-1)*dx; % spatial coordinates

% numerically integrate the pattern through z
pInt = p;%%%p is already sum across z dimension, then????%%%
%pInt = pInt(1:322,:);
%p = p(:,:,97)/Nx;

% Overal loop (more general than below since any background image can be used):
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For each spatial location:
%   Calculate displacements
%   For each line position/output image:
%       Add displaced pixels into image
%% 
image=rgb2gray(IMG_0064);
image=double(image);
step = round(size(image,2) / 1920);
imagesm=image(2100-160*step:step:2100+160*step,3300-221*step:step:3300+221*step);
%imagesm=imagesm/max(imagesm(:));
imagesm = imagesm / max(imagesm(:));
[imgs, d_x, d_y]= calcImgMC(pInt,Zd/dx*(k*dn0dp)/n0, imagesm);
%% 
% build a dictionary of signals at different amplitudes
dp = 10; % Pa, step size of dictionary
[D,plist] = build_dictionary(max(abs(pInt(:))),dp,Zd/dx*(k*dn0dp)/n0,omega,t);

% reconstruct the original pressure field
pRec = recon_p(imgs,D,plist);

% duplicate the reconstructed pressure field to all projection angles,
% and reconstruct the original field. 
maxAngle = 270; % degrees, maximum angle
dAngle = 1; % degrees
angles = 0:dAngle:maxAngle;
nAngles = length(angles);
pRecTomo = zeros(Nx,Ny,Ny);
for ii = 1:Nx
    pRecTomo(ii,:,:) = 1/dz*iradon(repmat(squeeze(pRec(ii,:)).',[1 nAngles]),angles,'v5cubic','Ram-Lak',1,Ny); %1/dz*iradon(repmat(sqz(pRec(ii,:)).',[1 nAngles]),angles,'v5cubic','Ram-Lak',1,Ny);
end
