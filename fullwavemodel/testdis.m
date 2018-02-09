apaz = apaz_sv(:,:,2:end);
z_sv = cumsum(zvec);
z_sv = z_sv(abs(z_sv-foc) < 0.25*foc);
num = find(abs(z_sv-foc) == min(abs(z_sv-foc)));
tTarg = (round((0.5*foc)/c0/dT):round((2*foc)/c0/dT)-1)*dT;

%% shift in z
apaz = reshape(apaz,[587,1,2881,116]);
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
% %%
n0 = 1.3325; % base index of refraction
% m, distance from center of US focus to background
Zd = 17.061e-2/2;
dn0dp = (1.3475-1.3325)/((1100-1)*100000); % dindex/dPascal, from Waxler paper, read from Figure 2 at 24.8 deg CSS
% %%
P = [];
tt = 1911:1953;
Pori = apaz_shifted(:,:,tt,:);
% P = Pori;
for ii = 1:size(Pori,4)
    P(:,:,:,ii) = permute(interp1(tTarg(tt),permute((Pori(:,:,:,ii)),[3 2 1]),tTarg(tt(1)):dT/4:tTarg(tt(end)),'PCHIP',0),[3 2 1]);%   apaz_shifted(:,:,ii) = interp1((t+tShift(ii))',squeeze(apaz_sv(:,:,ii))',tTarg','spline');%,0);
end
[x,y] = findpeaks(squeeze(P(round(size(P,1)/2),1,:,num)));
P = P(:,:,y(1):y(end),:);

%% interpolate 2D to 3D profile using interp1
cen = round([size(P,1),size(P,1)]/2);
[y1,x1] = meshgrid(1:587,1:587);
d = (sqrt((x1-cen(1)).^2+(y1-cen(2)).^2))*dX;
r =(-floor(nX/2):floor(nX/2))*dX;
d(d>293) = -1;
Pb = P;
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

dgrid = 92.6/98*1e-3; %m
% step = 1;
z_svnew = [flip(z_sv(num):-dgrid:z_sv(1)),z_sv(num)+dgrid:dgrid:z_sv(end)];
dxx = interp1((-floor(nX/2)+1:floor(nX/2)-1)*dX,dx(2:end-1,1,:,:),[flip(0:-dgrid:(-floor(nX/2)+1)*dX),dgrid:dgrid:(floor(nX/2)-1)*dX],'pchip');
dzz = interp1((-floor(nX/2)+1:floor(nX/2)-1)*dX,dz(2:end-1,1,:,:),[flip(0:-dgrid:(-floor(nX/2)+1)*dX),dgrid:dgrid:(floor(nX/2)-1)*dX],'pchip');
dzz = permute(interp1(z_sv(2:end-1),permute(dzz(:,1,:,2:end-1),[4 3 2 1]),z_svnew(2:end-1),'pchip'),[4 3 2 1]);
dxx = permute(interp1(z_sv(2:end-1),permute(dxx(:,1,:,2:end-1),[4 3 2 1]),z_svnew(2:end-1),'pchip'),[4 3 2 1]);
% step = step-1;

%%
img = imread('bkg.png');
bkg = double(img(2:end,2:end,1));
bkgnew = double(bkg);
row = squeeze(bkg(1,:));
rownum1 = find(row == 0);
col = squeeze(bkg(:,1));
colnum1 = find(col == 0);
row = squeeze(bkg(round(mean(rownum1(1:2))),:));
rownum2 = find(row == 0);
col = squeeze(bkg(:,round(mean(colnum1(1:2)))));
colnum2 = find(col == 0);

xx = (sort([colnum1;colnum2]));
zz = (sort([rownum1,rownum2]));



np =  mean(xx(3)-xx(1)); 
locxd = squeeze(round(dxx/(dgrid/(np/2))));
loczd = squeeze(round(dzz/(dgrid/(np/2))));
nphalf = np/2-1;
% np = np+2;
% nhalf = floor(np/2);
% np = np+1;
nt = size(locxd,2);
numx = size(locxd,1);
numz = size(loczd,3);
% normxz = max(max(abs(locxd(:))),max(abs(loczd(:))));
% locxd = flip(squeeze(locxd/normxz*nphalf),3);
% loczd = flip(squeeze(loczd/normxz*nphalf),3);
cen = floor((3*np)/2);
for i = 2:numx-10
    for j = 11:numz-10
        histmp = 255*ones(3*np,3*np,nt)/nt;
        v = bkg(xx(i)-nphalf:xx(i)+nphalf,zz(j)-nphalf:zz(j)+nphalf);
%         display([num2str(i) ' ' num2str(j) ' ' num2str(v) ' ']); 
        for kk = 1:nt
            x = round(cen+1+locxd(i,kk,j)); z = round(1+cen+loczd(i,kk,j));
            histmp(x-nphalf:x+nphalf,z-nphalf:z+nphalf,kk) = v/nt;
        end
        histmp = sum(histmp,3);
        bkgnew(xx(i)-nphalf:xx(i)+nphalf,zz(j)-nphalf:zz(j)+nphalf) = histmp(cen+1-nphalf:cen+1+nphalf,cen+1-nphalf:cen+1+nphalf);
    end
end
imshow(flip(uint8(bkgnew),2))