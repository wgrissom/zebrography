function [apaz_sv,z_sv,t,dT,nX,nY] = launch_asr_adaptive(p0,num,a,f0,c0)
omega0=2*pi*f0;
lambda=c0/f0; % wavelength
dX=lambda/num; dY=lambda/num;% grid spacing
nX=round(2*a*1.5/dX); nY=nX; % grid size
if(mod(nX,2)-1) % keep the grid odd
    nX=nX+1;
end
if(mod(nY,2)-1)
    nY=nY+1;
end
%  dZ=dX*5;
dZ = dX/2; % WAG quarter z-steps
% dT=dX/30/c0; %2*pi/omega0/20;
dT=dX/20/c0; %2*pi/omega0/20;
k=omega0/c0;
beta=8; % nonlinear coefficient
rho0=1000; % equilibrium density
N=beta/2/c0^3/rho0;
%p0=1e6; % transmitter pressure
%p0=(1e6)/4; % transmitter pressure % WAG quarter pressure
xaxis=(0:nX-1)*dX;xaxis=xaxis-mean(xaxis);
yaxis=(0:nY-1)*dY;yaxis=yaxis-mean(yaxis);
%% Generate circular aperture 
ap=zeros(nX,nY);
cen=round(size(ap)/2);
[i,j] = meshgrid(1:size(ap,1),1:size(ap,2));
ap = sqrt((i - cen(1)).^2 + (j - cen(2)).^2) <= round(a/dX);
% for i=1:nX
%     for j=1:nY
%         r=sqrt((i-cen(1))^2+(j-cen(2))^2);
%         if(r<=round(a/dX))
%             ap(i,j)=1;
%         end
%     end
% end

%%% Generate initial conditions based on input coordinates %%%%%%
ncycles = 8; % number of cycles in pulse
dur = 2; % exponential drop-off of envelope
duration=ncycles*6*2*pi/omega0;
nT=round(duration/dT);
if(mod(nT,2)-1)
    nT=nT+1;
end

t = (0:nT-1)*dT-2*ncycles/omega0*2*pi;
nSamples = round(ncycles/f0/dT);
hannWind = blackman(2*round(nSamples/2000)).';
window = [hannWind(1:length(hannWind)/2) ones(1,nSamples-length(hannWind)) ...
    hannWind(length(hannWind)/2+1:end)];
icvec = sin(t*omega0)*p0;
tWindow = (0:length(window)-1)*dT;
window = interp1(tWindow - mean(tWindow),window,t,'spline',0);
icvec = icvec.*window;
% icvec = exp(-(1.05*t*omega0/(ncycles*pi)).^(2*dur)).*sin(t*omega0)*p0;
%icvec = 1 ./ (1 + exp((abs(t) - ncycles/f0) ./ (10*ncycles/f0/128))) .*sin(t*omega0)*p0;

foc=2*a;
%
  [i,j] = meshgrid(1:size(ap,1),1:size(ap,2));
ix = i*dX - size(ap,1)/2*dX;
jy = j*dY - size(ap,2)/2*dY;
tt = sqrt(ix.^2+jy.^2+foc^2)/c0-(foc/c0);
tArray = permute(repmat(t(:),[1 size(ap)]),[2 3 1]) - repmat(tt,[1 1 nT]);
 ttWindow = permute(repmat(window(:),[1 size(ap)]),[2 3 1]) - repmat(tt,[1 1 nT]);
icvec = sin(tArray*omega0)*p0;
icvec = icvec.*ttWindow;
%  icvec = exp(-(1.05*(tArray)*omega0/(ncycles*pi)).^(2*dur)).*sin((tArray)*omega0)*p0;
% icvec = 1 ./ (1 + exp((abs(tArray) - ncycles/f0) ./ (10*ncycles/f0/128))) .*sin((tArray)*omega0)*p0;
apa = icvec.*repmat(ap,[1 1 nT]);
apa=flipdim(apa,3);

%% PRECALCULATE MODIFIED ANGULAR SPECTRUM 
[HH] = precalculate_mas(nX,nY,nT,dX,dY,dZ,dT,c0);
%% ABSORBING BOUNDARY LAYER
[abl] = precalculate_abl(nX,nY,nT);
%% PRECALCULATE ABSORPTION/DISPERSION FILTER
alpha0=2.17e-3; % dB/MHz^2/cm for water, note attenuation law proportional to f^2 here
[afilt3d] = precalculate_ad(alpha0,nX,nY,nT,dZ,dT);

%% MARCH ALONG PROPAGATION AXIS %%
flag = 0;
apaz=apa; zvec=dZ; dZa=dZ; cc=1;
kkk =1;
apaz_sv = [];
z_sv = [];
while(sum(zvec)<1.5*foc)
    dZa=dZ;
    %N*dZ/dT*max(max(max(apaz)))
    if(N*dZa/dT*max(max(max(apaz)))>0.1) % WAG quarter z-steps
    %if(N*dZ/dT*max(max(max(apaz)))>0.1)
        %dZa=0.1*dT/max(max(max(apaz)))/N;
        dZa=0.1*dT/max(max(max(apaz)))/N; % WAG quarter z-steps
    end
    zvec(cc)=dZa;
    disp(['Propagation distance = ' num2str(sum(zvec)) ' m'])
    disp(['Current z-step = ' num2str(dZa) ' m'])
%     
    apaz=march_asr(apaz,dZa,dT,N,HH,abl,afilt3d); 
     if abs(sum(zvec(:))-foc) < 15*dZ %6
%     if abs(zvec(cc)-foc) < 4*dZ
        apaz_sv(:,:,:,kkk) = apaz(1:end,1:end,:); % the x/y/time profile at current z-coordinate
        z_sv(kkk) = sum(zvec); % the current z-coordinate
        kkk = kkk+1;
        flag = 1;
     elseif flag == 1
        break;
     end
%     if length(z_sv)>1 && abs(sum(zvec(:))-foc) > 5*dZ
%         break;
%     end
    %apaz_sv(:,:,end+1) = squeeze(apaz(round(nX/2),:,:));
    cc=cc+1;
end
clear apaz