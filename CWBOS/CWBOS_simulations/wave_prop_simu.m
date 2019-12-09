function [apaz_sv,dX,dY,dZ,nX,nY,z_sv] = wave_prop_simu(p0,a,f0,f_num,is_save)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% WRITTEN BY GIANMARCO PINTON 
% ON 2017-05-01
% NONLINEAR PROPAGATION IN HOMOGENEOUS MEDIUM
% MODIFIED ANGULAR SPECTRUM 
% RUSANOV
% FREQUENCY DOMAIN ATTENUATION AND DISPERSION
% ABSORBING BOUNDARY LAYERS
% ADAPTIVE STEP SIZE IN PROPAGATION 
% THIS VERSION DESIGNED FOR WILL GRISSOM AND HUIWEN LUO

%| Inputs:
%|  p0     Transmitter pressure; set to match optical hydrophone;
%|  f0     Center frequency of the transducer
%|  a      source radius of the transducer
%|  f_num  f-number of the transducer
%|  is_save save all the data in simulations in the current path
%|
%| Outputs:
%| apaz_sv   [Nx,Nt,Nz] simulated spatially and temporally-resolved pressure
%| dX,dY .   Grid spacing, which is set as 1/8 of wavelength in this function.
%| dZ        Wave propagation's step size 
%|
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
beta=3.5; % nonlinear coefficient in distilled water
rho0=1000; % equilibrium density
widePulse = true;
omega0=2*pi*f0;
c0=1500; % speed of sound
lambda=c0/f0; % wavelength
dX=lambda/8; dY=lambda/8;% grid spacing
nX=round(2*a*1.5/dX); nY=nX; % grid size
if(mod(nX,2)-1) % keep the grid odd
    nX=nX+1;
end
if(mod(nY,2)-1)
    nY=nY+1;
end
dZ=dX*2;%  z-steps
dT=dX/5/c0; 
k=omega0/c0;
N=beta/2/c0^3/rho0;
xaxis=(0:nX-1)*dX;xaxis=xaxis-mean(xaxis);
yaxis=(0:nY-1)*dY;yaxis=yaxis-mean(yaxis);

%%% Generate circular aperture
ap=zeros(nX,nY);
cen=round(size(ap)/2);
[i,j] = meshgrid(1:size(ap,1),1:size(ap,2));
ap = sqrt((i - cen(1)).^2 + (j - cen(2)).^2) <= round(a/dX);


%%% Generate initial conditions based on input coordinates %%%%%%
if widePulse
    ncycles = 12; % number of cycles in pulse (originally 2)
else
    ncycles = 2;
end
dur = 2; % exponential drop-off of envelope
duration=ncycles*6*2*pi/omega0;
nT=round(duration/dT);
if(mod(nT,2)-1)
    nT=nT+1;
end


t = (0:nT-1)*dT-2*ncycles/omega0*2*pi;

if widePulse
    nSamples = round(ncycles/f0/dT);
    hannWind = blackman(1.5*round(nSamples/10)).'; %10000
    window = [hannWind(1:length(hannWind)/2) ones(1,nSamples-length(hannWind)) ...
        hannWind(length(hannWind)/2+1:end)];
    icvec = sin(t*omega0)*p0;
    tWindow = (0:length(window)-1)*dT;
    window = interp1(tWindow - mean(tWindow),window,t,'spline',0);
    icvec = icvec.*window;
else
    icvec = exp(-(1.05*t*omega0/(ncycles*pi)).^(2*dur)).*sin(t*omega0)*p0;
end

foc=f_num*a;
%
[i,j] = meshgrid(1:size(ap,1),1:size(ap,2));
ix = i*dX - size(ap,1)/2*dX;
jy = j*dY - size(ap,2)/2*dY;
tt = sqrt(ix.^2+jy.^2+foc^2)/c0-(foc/c0);
tArray = permute(repmat(t(:),[1 size(ap)]),[2 3 1]) - repmat(tt,[1 1 nT]);
icvec = exp(-(1.05*(tArray)*omega0/(ncycles*pi)).^(2*dur)).*sin((tArray)*omega0)*p0;
apa = icvec.*repmat(ap,[1 1 nT]);
clear icvec tArray
apa=flip(apa,3);
%% PRECALCULATE MODIFIED ANGULAR SPECTRUM
[HH] = precalculate_mas(nX,nY,nT,dX,dY,dZ,dT,c0);
%% ABSORBING BOUNDARY LAYER
[abl] = precalculate_abl(nX,nY,nT);
%% PRECALCULATE ABSORPTION/DISPERSION FILTER
alpha0=2.17e-3; % dB/MHz^2/cm for water, note attenuation law proportional to f^2 here
[afilt3d] = precalculate_ad(alpha0,nX,nY,nT,dZ,dT);

%% MARCH ALONG PROPAGATION AXIS %%
apaz=apa; zvec=dZ; cc=1;clear apa
apaz_sv = [];
z_sv = [];
kk = 1;
while sum(zvec) < 1.5*foc
    dZa=dZ;
    if(N*dZ/dT*max(max(max(apaz)))>0.1) % WAG quarter z-steps
        dZa=0.1*dT/max(max(max(apaz)))/N; % WAG quarter z-steps
    end
    zvec(cc)=dZa;
    disp(['Propagation distance = ' num2str(sum(zvec)) ' m'])
    disp(['Current z-step = ' num2str(dZa) ' m'])
    apaz=march_asr(apaz,dZa,dT,N,HH,abl,afilt3d);
    if abs(sum(zvec)-foc) < 0.5*foc
        apaz_sv(:,:,kk) = squeeze(apaz(round(nX/2),:,:));
        z_sv(kk) = sum(zvec(:));
        kk = kk + 1;
    end
    cc=cc+1;
end
clear apaz HH abl afilt3d
apaz_sv = single(apaz_sv);
if is_save
    save(['p0_',num2str(p0),'_',num2str(f_num),'.mat'], '-v7.3');
end
