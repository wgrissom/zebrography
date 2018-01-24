function [hisdic,projdic,P] = generate_dict(apaz_shifted,z_sv,dX,dT,tTarg,num)
%num  slice number of the focus

%% build the dictionary (historgrams -- projected waveforms)

%%
n0 = 1.3325; % base index of refraction
%  Zd = dY*nY/2; % m, distance from center of US focus to background
Zd = 0.05;
dn0dp = (1.3475-1.3325)/((1100-1)*100000); % dindex/dPascal, from Waxler paper, read from Figure 2 at 24.8 deg CSS
%%
P = [];
% tt = 715:745;
% tt = 1681:1723;
tt = 635:653;
 Pori = apaz_shifted(:,:,tt,:);
% P = Pori;
for ii = 1:size(Pori,4)
    P(:,:,:,ii) = permute(interp1(tTarg(tt),permute(squeeze(Pori(:,:,:,ii)),[3 2 1]),tTarg(tt(1)):dT/10:tTarg(tt(end)),'PCHIP',0),[3 2 1]);%   apaz_shifted(:,:,ii) = interp1((t+tShift(ii))',squeeze(apaz_sv(:,:,ii))',tTarg','spline');%,0);
end
[x,y] = findpeaks(squeeze(P(round(size(P,1)/2),round(size(P,1)/2),:,num)));
P = P(:,:,y(1):y(end)-1,:);
% ppp = permute(P,[3,1,2,4]);
% ppp = downsample(ppp,9);
% P = permute(ppp,[2,3,1,4]);
% P = Pori;
%% project pressure waveforms
% calculate the projections using Gmri
angles = 0:179;
kk = -1:2/(size(P,1)-1):1;
kk = [kron(cosd(angles(:)),kk.') kron(sind(angles(:)),kk.')];
G = Gmri(kk,true(size(P,1)),'fov',size(P,1)/2);
proj = zeros([size(P,1),length(angles),size(P,3),size(P,4)]);
for i = 1:length(z_sv)
    for j = 1:size(P,3)
        I = reshape(G * (squeeze(P(:,:,j,i))),[size(P,1),size(angles)]);
        proj(:,:,j,i) = real(ifftshift(ifft(ifftshift(I,1),[],1)));
    end
end

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

% dx = squeeze(Zd*1/n0*dn0dp*Ix);
%% generate 2D histogram at each x-z location
locz = zeros(size(dz));
flag1 = dz>0;
flag2 = dz<0;
for i = 2:size(proj,4)-1
    tmpdz = dz(:,:,:,i);
    tmplocz = locz(:,:,:,i);
    tmplocz(flag1(:,:,:,i)) = tmpdz(flag1(:,:,:,i))/(z_sv(i+1)-z_sv(i));
    tmplocz(flag2(:,:,:,i)) = tmpdz(flag2(:,:,:,i))/(z_sv(i)-z_sv(i-1));
%     locz(dz>0) = dz(dz>0)/(z_sv(i+1)-z_sv(i));
%     locz(dz<0) = dz(dz<0)/(z_sv(i)-z_sv(i-1));
    locz(:,:,:,i) = tmplocz;
end
locx = dx/dX;

locx = permute(locx,[2 3 1 4]);
locx = locx(:,:,:);
locz = permute(locz,[2 3 1 4]);
locz = locz(:,:,:);

%%
xxx = round(max(abs(locz(:))));
yyy = round(max(abs(locx(:))));
xxx = max(xxx,yyy);
yyy = xxx;
hissize = 15;
his = zeros(hissize,hissize,size(locx,1), size(locx,3));
histmp = zeros(2*yyy+1,2*xxx+1);
for i = 1:size(locx,1)
    for j = 1:size(locx,3)
        for kk = 1:size(locx,2)
            x = round(yyy+1+locx(i,kk,j)); z = round(xxx+1+locz(i,kk,j));
            histmp(x,z)=histmp(x,z)+1;
        end
        if length(find(histmp ~=0)) == 1
            his(:,:,i,j) = zeros(hissize,hissize);
            his(round(hissize/2),round(hissize/2),i,j) = 1;
        else
            tmp = imresize(histmp/size(locz,2),[hissize,hissize]);
            tmp(tmp<0) = 0;
            his(:,:,i,j) = tmp/sum(tmp(:));
        end
        histmp = zeros(2*yyy+1,2*xxx+1);
    end
end

%%
% ghis = zeros(size(his,3)*hissize,size(his,4)*hissize);
% for i = 1:size(locx,1)
%     for j = 1:size(locx,3)
%         ghis(hissize*(i-1)+1:hissize*i,hissize*(j-1)+1:hissize*j) = his(:,:,i,j);
%     end
% end

hisdic = his;
projdic = permute(proj,[2 3 1 4]);
projdic = projdic(:,:,:);
% hisdic = permute(hisdic,[3 1 2]);
