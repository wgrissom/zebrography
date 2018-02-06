function [locx,locz,P] = generate_simu(apaz_shifted,z_sv,dX,dT,tTarg,num)
%num  slice number of the focus
%bkgdot original background pattern  
%dgrid grid spacing in the real experiment.
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
% tt = 635:653;
% tt = 387:397;  %newdata
% tt = 443:457; %latestdata
% tt = 290:302;
% tt = 1751:1793; %013118
% tt = 1912:1954;
tt = 1911:1953;
Pori = apaz_shifted(:,:,tt,:);
% P = Pori;
for ii = 1:size(Pori,4)
    P(:,:,:,ii) = permute(interp1(tTarg(tt),permute((Pori(:,:,:,ii)),[3 2 1]),tTarg(tt(1)):dT/50:tTarg(tt(end)),'PCHIP',0),[3 2 1]);%   apaz_shifted(:,:,ii) = interp1((t+tShift(ii))',squeeze(apaz_sv(:,:,ii))',tTarg','spline');%,0);
end
% [x,y] = findpeaks(squeeze(P(round(size(P,1)/2),round(size(P,1)/2),:,num)));
[x,y] = findpeaks(squeeze(P(round(size(P,1)/2),1,:,num)));
P = P(:,:,y(1):y(end),:);
% ppp = permute(P,[3,1,2,4]);
% ppp = downsample(ppp,9);
% P = permute(ppp,[2,3,1,4]);
% P = Pori;
%% project pressure waveforms
% calculate the projections using Gmri
% angles = 0;  %% changed it to one angle 
% kk = -1:2/(size(P,1)-1):1;
% kk = [kron(cosd(angles(:)),kk.') kron(sind(angles(:)),kk.')];
% G = Gmri(kk,true(size(P,1)),'fov',size(P,1)/2);
% proj = zeros([size(P,1),length(angles),size(P,3),size(P,4)]);
% for i = 1:length(z_sv)
%     for j = 1:size(P,3)
%         I = reshape(G * (squeeze(P(:,:,j,i))),[size(P,1),size(angles)]);
%         proj(:,:,j,i) = real(ifftshift(ifft(ifftshift(I,1),[],1)));
%     end
% end
proj = P;
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

% dx = squeeze(Zd*1/n0*dn0dp*Ix);
%% generate 2D histogram at each x-z location
% locz = zeros(size(dz));
% flag1 = dz>0;
% flag2 = dz<0;
dZ = z_sv(2)-z_sv(1);
% for i = 2:size(proj,4)-1
%     tmpdz = dz(:,:,:,i);
%     tmplocz = locz(:,:,:,i);
%     tmplocz(flag1(:,:,:,i)) = tmpdz(flag1(:,:,:,i))/((z_sv(i+1)-z_sv(i)));
%     tmplocz(flag2(:,:,:,i)) = tmpdz(flag2(:,:,:,i))/((z_sv(i)-z_sv(i-1)));
% %     locz(dz>0) = dz(dz>0)/(z_sv(i+1)-z_sv(i));
% %     locz(dz<0) = dz(dz<0)/(z_sv(i)-z_sv(i-1));
%     locz(:,:,:,i) = tmplocz;
% end
nX = size(apaz_shifted,1);
dgridb = 92.6/98*1e-3; %m
step = 10;
dgrid = dgridb/step;
z_svnew = [flip(z_sv(40):-dgrid:z_sv(1)),z_sv(40)+dgrid:dgrid:z_sv(end)];
dxx = interp1((-floor(nX/2)+1:floor(nX/2)-1)*dX,dx(2:end-1,1,:,:),[flip(0:-dgrid:(-floor(nX/2)+1)*dX),dgrid:dgrid:(floor(nX/2)-1)*dX],'pchip');
dzz = interp1((-floor(nX/2)+1:floor(nX/2)-1)*dX,dz(2:end-1,1,:,:),[flip(0:-dgrid:(-floor(nX/2)+1)*dX),dgrid:dgrid:(floor(nX/2)-1)*dX],'pchip');
dzz = permute(interp1(z_sv(2:end-1),permute(dzz(:,1,:,2:end-1),[4 3 2 1]),z_svnew(2:end-1),'pchip'),[4 3 2 1]);
dxx = permute(interp1(z_sv(2:end-1),permute(dxx(:,1,:,2:end-1),[4 3 2 1]),z_svnew(2:end-1),'pchip'),[4 3 2 1]);
locz = dzz/dgrid;
locx = dxx/dgrid;
locx = locx([flip(round(size(dxx,1)/2):-step:1),round(size(dxx,1)/2)+step:step:end],:,:,[flip(round(length(z_svnew)/2):-step:1),round(length(z_svnew)/2)+step:step:end])/step;
locz = locz([flip(round(size(dxx,1)/2):-step:1),round(size(dxx,1)/2)+step:step:end],:,:,[flip(round(length(z_svnew)/2):-step:1),round(length(z_svnew)/2)+step:step:end])/step;
% locx = dx/dX;
% locz = dz/dZ;
locx = squeeze(permute(locx,[2 3 1 4]));
% locx = locx(:,:,:);
locz = squeeze(permute(locz,[2 3 1 4]));
% locz = locz(:,:,:);

% %%
% 
% numz = size(locx,4)-2;
% numx = size(locx,3)-2;
% np = 11;
% nphalf = floor(np/2);
% img = zeros(numx*np,numz*np);
% % img(np*((1:numx))-nphalf,np*((1:numz))-nphalf) = 1;
% img(np*((1:2:numx))-nphalf,np*((1:2:numz))-nphalf) = 1;
% img(np*((2:2:numx))-nphalf,np*((2:2:numz))-nphalf) = 1;
% 
% %%
% normxz = max(max(abs(locx(:))),max(abs(locz(:))));
% locx = squeeze(locx/normxz*nphalf);
% locz = squeeze(locz/normxz*nphalf);
% 
% for i = 1:numx
%     for j = 1:numz
%         histmp = zeros(np,np);
%         v = img(np*i-nphalf,np*j-nphalf);
%         display([num2str(i) ' ' num2str(j) ' ' num2str(v) ' ']); 
%         for kk = 1:nt
%             x = round(1+nphalf+locx(kk,i,j)); z = round(1+nphalf+locz(kk,i,j));
%             histmp(x,z) = histmp(x,z)+v/nt;
%         end
%         img((i*np-np+1):(i*np),(j*np-np+1):(j*np)) = histmp;
%     end
% end
% 
% % for i = 2:2:numz
% %     for j = 2:2:numx
% %         histmp = zeros(np,np);
% %         for kk = 1:nt
% %             x = round(nphalf+1+locx(i,kk,j)); z = round(nphalf+1+locz(i,kk,j));
% %             histmp(x,z) = histmp(x,z)+1;
% %         end
% %         img(i-nphalf:i+nphalf,j-nphalf:j+nphalf) = histmp/sum(histmp(:));
% %     end
% % end
% 
% % for i = 1
% %     tic
% %     for j = 1:size(locx,3)
% %         for kk = 1:size(locx,2)
% %             x = round(yyy+1+locx(i,kk,j)); z = round(xxx+1+locz(i,kk,j));
% %             histmp(x,z)=histmp(x,z)+1;
% %         end
% % %         x = squeeze(round(yyy+1+locx(i,:,j)));
% % %         z = squeeze(round(xxx+1+locz(i,:,j)));
% %         tmp = zeros(2*yyy+1,2*xxx+1,size(locx,2));
% %         tmp(x(:),z(:),1:size(locx,2)) = 1;
% %         histmp = sum(tmp,3);
%         if length(find(histmp ~=0)) == 1
%             his(:,:,i,j) = zeros(hissize,hissize);
%             his(round(hissize/2),round(hissize/2),i,j) = 1;
%         else
%             tmp = histmp/size(locz,2);
%             tmp(tmp<0) = 0;
%             his(:,:,i,j) = tmp/sum(tmp(:));
%         end
%         histmp = zeros(2*yyy+1,2*xxx+1);
%     end
%     toc
% end

%%his = squeeze(his(:,:,1,:));
% his = reshape(his,15,15,75,60);
% size(his)
% ghis = zeros(size(his,3)*hissize,size(his,4)*hissize);
% for i = 1:75
% for j = 1:60
% ghis(hissize*(i-1)+1:hissize*i,hissize*(j-1)+1:hissize*j) = his(:,:,i,j);
% end
% end

%%
% ghis = zeros(size(his,3)*hissize,size(his,4)*hissize);
% for i = 1:size(locx,1)
%     for j = 1:size(locx,3)
%         ghis(hissize*(i-1)+1:hissize*i,hissize*(j-1)+1:hissize*j) = his(:,:,i,j);
%     end
% end
% 
% hisdic = his;
% projdic = permute(proj,[2 3 1 4]);
% projdic = projdic(:,:,:);
% hisdic = permute(hisdic,[3 1 2]);
