%generate simuated image for the real image 
%locx and locz are displacements from genrate_simu

locxd = locx;%(:,[flip(294:-6:1),294+6:6:end],[flip(68:-3:1),68+3:3:end])/(dgrid/dX);
loczd = locz;%(:,[flip(294:-6:1),294+6:6:end],[flip(68:-3:1),68+3:3:end])/(dgrid/dZ);
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
zz = flip(sort([rownum1,rownum2]));

np =  mean(xx(2:end)-xx(1:end-1))*3                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           ;
nphalf = np/2;
np = np+1;
nt = size(locx,1);
numx = size(locxd,2);
numz = size(loczd,3);
normxz = max(max(abs(locxd(:))),max(abs(loczd(:))));
locxd = flip(squeeze(locxd/normxz*(nphalf)),3);
loczd = flip(squeeze(loczd/normxz*(nphalf)),3);

for i = 10:numx-10
    for j = 10:numz-10
        histmp = 255*ones(2*np-1,2*np-1,nt)/nt;
        v = bkg(xx(i)-nphalf:xx(i)+nphalf,zz(j+14)-nphalf:zz(j+14)+nphalf);
%         display([num2str(i) ' ' num2str(j) ' ' num2str(v) ' ']); 
        for kk = 1:nt
            x = round(2*nphalf+1+locxd(kk,i,j)); z = round(1+2*nphalf+loczd(kk,i,j));
            histmp(x-nphalf:x+nphalf,z-nphalf:z+nphalf,kk) = v/nt;
        end
        histmp = sum(histmp,3);
        bkgnew(xx(i)-nphalf+1:xx(i)+nphalf-1,zz(j+14)-nphalf+1:zz(j+14)+nphalf-1) = histmp(nphalf+2:2*np-2-nphalf,nphalf+2:2*np-2-nphalf);
    end
end
imshow(uint8(bkgnew),[]);
