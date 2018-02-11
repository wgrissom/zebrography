% apaz = apaz_sv(:,:,2:end);
z_sv = cumsum(zvec);
ind = find(abs(z_sv-foc) < 0.5*foc);
apaz = apaz_sv(:,:,ind);
z_sv = z_sv(ind);
num = find(abs(z_sv-foc) == min(abs(z_sv-foc)));
tTarg = (round((0.5*foc)/c0/dT):round((2*foc)/c0/dT)-1)*dT;

%% shift in z
apaz = reshape(apaz,[587,1,2881,size(apaz,3)]);
zPos = z_sv;
tShift = zPos/c0; % time shifts (position/(m/s) = s)
apaz(isnan(apaz)) = 0;

%%
apaz_shifted = [];
for ii = 1:size(apaz,4)
     apaz_shifted(:,:,:,ii) = permute(interp1((t+tShift(ii))',permute((apaz(:,:,:,ii)),[3 2 1 ]),tTarg','spline',0),[3 2 1]);
%     apaz_shifted(:,:,ii) = interp1((t+tShift(ii))',squeeze(apaz_sv(:,:,ii))',tTarg','spline');%,0);
    
end
% 
clear apaz_sv
% %%
n0 = 1.3325; % base index of refraction
% m, distance from center of US focus to background
Zd = 17.061e-2/2;
dn0dp = (1.3475-1.3325)/((1100-1)*100000); % dindex/dPascal, from Waxler paper, read from Figure 2 at 24.8 deg CSS
% %%
P = [];
tt = 1896:1938;
Pori = apaz_shifted(:,:,tt,:);
% P = Pori;
for ii = 1:size(Pori,4)
    P(:,:,:,ii) = permute(interp1(tTarg(tt),permute((Pori(:,:,:,ii)),[3 2 1]),tTarg(tt(1)):dT/4:tTarg(tt(end)),'PCHIP',0),[3 2 1]);%   apaz_shifted(:,:,ii) = interp1((t+tShift(ii))',squeeze(apaz_sv(:,:,ii))',tTarg','spline');%,0);
end
[x,y] = findpeaks(squeeze(P(round(size(P,1)/2),1,:,num)));
P = P(:,:,y(1):y(end),:);
clear apaz_shifted;
clear Pori
%% interpolate 2D to 3D profile using interp1
cen = round([size(P,1),size(P,1)]/2);
[y1,x1] = meshgrid(1:587,1:587);
d = (sqrt((x1-cen(1)).^2+(y1-cen(2)).^2))*dX;
r =(-floor(nX/2):floor(nX/2))*dX;
d(d>293) = -1;
Pb = P;
clear P;
Pnew = zeros(size(Pb,1),size(Pb,1),size(Pb,3),size(Pb,4));
for ii = 1:size(Pnew,3)
    for jj = 1:size(Pnew,4)
    Pnew(:,:,ii,jj) = interp1(r,squeeze(Pb(:,1,ii,jj)),d,'spline',0);
    end
end
%%
n0 = 1.3325; % base index of refraction
% m, distance from center of US focus to background
Zd = 17.061e-2/2;
dn0dp = (1.3475-1.3325)/((1100-1)*100000); % dindex/dPascal, from Waxler paper, read from Figure 2 at 24.8 deg CSS
proj = sum(Pnew,2)*dX;
nt = size(proj,3);
dz = zeros(size(proj));
% loop through z-positions and calculated centered differences
for i = 2:length(z_sv)-1
    dz(:,:,:,i) = Zd*1/n0*dn0dp* (proj(:,:,:,i+1)-proj(:,:,:,i-1))/((z_sv(i+1)-z_sv(i-1)));
    
end

%% caculate displacement waveforms in x-z plane
% dz = Zd*1/n0*dn0dp*Iz;

dx = zeros(size(proj));
% loop through x-positions and calculated centered differences
for i = 2:size(proj,1)-1
    dx(i,:,:,:) = Zd*1/n0*dn0dp*(proj(i+1,:,:,:)-proj(i-1,:,:,:))/(2*dX);
end

%% for real img
ds = 3e-3/16/8;

locx = (squeeze(round(dx/(ds))));
locz = (squeeze(round(dz/(ds))));

locxd = flip(locx,3);   %displacements for locations
loczd = -flip(locz,3);
% locxd = 1.2*locx;
% loczd = 1.2*locz;
nz = size(loczd,3);
nt = size(loczd,2);
xloc = kron(ones(5,1),1:587);
zloc = kron(ones(5,1),1:nz);
xloc = xloc(:);
zloc = zloc(:);
[zo,xo] = meshgrid(1:nz,1:587);
[zz,xx] = meshgrid(zloc,xloc);
% [xo,zo] = meshgrid(1:587,1:nz);
% [xx,zz] = meshgrid(xloc,zloc);
for ii = 1:nt
    locxx(:,ii,:) = interp2(zo,xo,squeeze(locxd(:,ii,:)),zz,xx);
    loczz(:,ii,:) = interp2(zo,xo,squeeze(loczd(:,ii,:)),zz,xx);
end
% locxx = flip(locxx,3);
[zzz,xxx] = meshgrid(1:length(zloc),1:length(xloc));
xdisloc = round(repmat(reshape(xxx,size(xxx,1),1,size(xxx,2)),1,nt,1)+locxx);
zdisloc = round(repmat(reshape(zzz,size(zzz,1),1,size(zzz,2)),1,nt,1)+loczz);
%caculate coordinates for each pixel at each time point

%%

img = imread('./2018_02_09/IMG_0023.CR2');
histmp = zeros(size(xdisloc,1),size(xdisloc,3),3);
partimg = double(img(1+1026:length(xloc)+1026,1+2386+43:length(zloc)+2386+43,:));
for ii = 1:nt
    ind = sub2ind([size(xdisloc,1),size(zdisloc,3)],squeeze(xdisloc(:,ii,:)),squeeze(zdisloc(:,ii,:)));
    for jj = 1:3
        tmp = partimg(:,:,jj);
        histmp(:,:,jj)= histmp(:,:,jj) + double(tmp(ind));
    end
end
bkg = img;
bkg(1+1026:length(xloc)+1026,1+2386+43:length(zloc)+2386+43,:) = round(histmp/nt);

figure; imshow(uint8(round(bkg(:,1+2386+43:length(zloc)+2386+43,:))));
saveas(gcf,'./simuing/simuimg_0023.jpg');   