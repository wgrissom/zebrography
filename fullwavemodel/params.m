p0 = 300/160.94*1e6;
a=63.2/2/1e3;%%%  % source radius
f0=1.1e6; % center frequency of pulse
omega0=2*pi*f0;
c0=1500; % speed of sound
foc = 2*a;
k = 2*pi*f0/c0;
omega0=2*pi*f0;
lambda=c0/f0; % wavelength
num = 2;
dX=lambda/num; dY=lambda/num;% grid spacing
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
dT=dX/5/c0; %2*pi/omega0/20;
k=omega0/c0;
cd ../kzk_toHuiwen/irt;
setup;
cd ..;
cd ../fullwavemodel;
