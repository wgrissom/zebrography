
%% pressure pre-processing 
apaz = apaz_sv; %dimension:[nX,nT,nZ]
nX = size(apaz,1);
nT = size(apaz_sv,2);
apaz = reshape(apaz,[nX,1,nT,size(apaz,3)]);
apaz(isnan(apaz)) = 0;
clear apaz_sv
apaz = squeeze(apaz);
apaz(:,size(apaz,2)+1:4500,:) = 0; % zero padding for shifting in z

sig = squeeze(max(max(apaz,[],2),[],1));
nn = find(sig == max(sig));
display(nn); % find the slice where the focus is.(Peak postive 
z_svnew = [flip(z_sv(nn):-dZ:z_sv(1)),z_sv(nn)+dZ:dZ:z_sv(end)];% 
apaznew = permute(interp1(z_sv,permute(apaz,[3 2 1]),z_svnew,'PCHIP',0),[3 2 1]); 
%Interpolate pressure field along z-dimension
%shift in z-dimension  
for ii = 1:size(apaznew,3)
    apaz_shifted(:,:,ii) = circshift(squeeze(apaznew(:,:,ii)),[0 10*ii]); 
    %dZ = 2*dX
    %dT = dX/5/c0
    %dZ = 10*(dT*c0)
end

% find an ultrasound cycle 
sig = squeeze(max(max(apaz_shifted,[],2),[],1));
nn = find(sig == max(sig));
tmp = squeeze(apaz_shifted(round(nX/2),:,nn));
tt = find(tmp == max(tmp(:))):find(tmp == max(tmp(:)))+40;
% tt = 2739:2779
P = apaz_shifted(:,tt,:);

% interpolate 2D to 3D profile using interp1
P = squeeze(P);
P = reshape(squeeze(P),size(P,1),1,size(P,2),size(P,3));%dim:[nX,1,nT,nZ]
cen = round([size(P,1),size(P,1)]/2);
[y1,x1] = meshgrid(1:nX,1:nX);
d = (sqrt((x1-cen(1)).^2+(y1-cen(2)).^2));
r =(-floor(nX/2):floor(nX/2));
d(d>293) = -inf;
Pb = P;
Pnew = zeros(size(Pb,1),size(Pb,1),size(Pb,3),size(Pb,4));% dim:[nX,nX,nT,nZ]
% Pnew = zeros(size(Pb,1),size(Pb,1),size(Pb,4),size(Pb,3));
for ii = 1:size(Pb,3) % time dimension
    for jj = 1:size(Pb,4) % z dimension
%       Pnew(:,:,jj,ii) = interp1(r,squeeze(Pb(:,1,ii,jj)),d,'spline',0);
    Pnew(:,:,ii,jj) = interp1(r,squeeze(Pb(:,1,ii,jj)),d,'spline',0);
    end
end

%project pressure along line-of-sight dimension 
proj = sum(Pnew,2)*dX;% Pnew dim:[nX,nX,nT,nZ] proj: [nX,1,nT,nZ]

%% central finite differencing to get displacements 
n = 4; %accuracy
A = zeros(n+1,n+1);
for ii = 1:n+1
    A(ii,:) = (-n/2:n/2).^(ii-1);
end
c = A\[0;1;zeros(n-1,1)];
n0 = 1.3325; % base index of refraction
% m, distance from center of US focus to background
Zd = 17.061e-2/2;
dn0dp = (1.3475-1.3325)/((1100-1)*100000); % dindex/dPascal, from Waxler paper, read from Figure 2 at 24.8 deg CSS
 % integrate through line of sight
dz = zeros(size(proj));
% loop through z-positions and calculated centered differences
nt = size(proj,3);
for kk = 1:nt
    for ii = 1+n/2:length(z_svnew)-1-n/2 % z dim
     dz(:,1,kk,ii) = Zd*1/n0*dn0dp*squeeze(proj(:,1,kk,ii-n/2:ii+n/2))*c/dZ; %proj(:,z-2:z+2)*c
    end
end
dx = zeros(size(proj));
% loop through x-positions and calculated centered differences
for kk = 1:nt
    for ii = 1+n/2:size(proj,1)-1-n/2
        dx(ii,1,kk,:) = Zd*1/n0*dn0dp*c'*squeeze(proj(ii-n/2:ii+n/2,1,kk,:))/dX; %c' * proj(:,x-2:x+2)
    end
end

%% apply displacements to unblurred photo to get simulated histograms

ds = 3e-3/16/7; % the distance of each pixel in the real photo 
nn1 = round(dX/ds);
nn2 = round(dZ/ds); 
locx = (squeeze(round(dx/(ds))));
locz = (squeeze(round(dz/(ds))));
locxd = locx;   %displacements for locations
loczd = locz;
offsetx = 212;
offsetz = 114;
nt = size(locx,2);
nx = size(locx,1);
nz = size(locx,3);
img = rgb2gray(imread(['./2018_02_09/IMG_0019.CR2'])); % convert RGB image to gray image
%take the samller image field for the simulated pressure
partimg = double(img(1+offsetx:nx*nn1+offsetx,1+offsetz:nz*nn2+offsetz,:)); % just select the good part used for simulation from the real photo 
block = partimg(1820-50:1820+50,1991-50:1991+50,:);
% imshow(uint8(block));
pin = block(51-5:51+5,51-5:51+5,:);  % find pin 
N = round(size(block,1)); 
nn = 5;
hissize = N;
bkg = zeros(N,N,1);
[x,y] = meshgrid(1:43);
[xq,yq] = meshgrid(linspace(1,43,N));
for ii = 1:1
    bkg(:,:,ii) = interp2(x,y,block(1:43,1:43),xq,yq,'spline');
end
bkg = imgaussfilt(bkg,4);
bkg = ones(size(bkg))*mean(bkg(:)); %smooth background image 
block = bkg;
block(51-5:51+5,51-5:51+5) = pin;
imshow(uint8(pin));

%% histogram for pressure of each point

his = zeros(hissize,hissize,1,nx,nz);
histmp = repmat(bkg,1,1,1,nt);
for ii = 1:nx%294-60:294+60
    for jj =1:nz%112-60:112+60
        for kk = 1:nt
            x = round(round(N/2)+locxd(ii,kk,jj)); z = round(round(N/2)+loczd(ii,kk,jj));
            histmp(x-nn:x+nn,z-nn:z+nn,:,kk)= pin;
        end 
%         figure;
%         imshow(uint8(sum(histmp,4)/size(locx,2)));
        his(:,:,:,ii,jj) = sum((histmp),4)/nt;
        histmp = repmat(bkg,1,1,1,nt);
    end
end
