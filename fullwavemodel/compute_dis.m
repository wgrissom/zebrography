function [dz,dx,P] = compute_dis(apaz_shifted,z_sv,dX,dT,tTarg,num)
%num  slice number of the focus
%bkgdot original background pattern  
%dgrid grid spacing in the real experiment.

n0 = 1.3325; % base index of refraction
% m, distance from center of US focus to background
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
%tt = 1911:1953;
tt = 1889:1931;
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
