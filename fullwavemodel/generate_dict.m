function [hisdic,projdic] = generate_dict(apaz_shifted,z_sv,dX,dT,tTarg)

%% build the dictionary (historgrams -- projected waveforms)

%%
n0 = 1.3325; % base index of refraction
%  Zd = dY*nY/2; % m, distance from center of US focus to background
 Zd = 0.05;
dn0dp = (1.3475-1.3325)/((1100-1)*100000); % dindex/dPascal, from Waxler paper, read from Figure 2 at 24.8 deg CSS
%%
% Pori = apaz_shifted(:,:,:,188-10:188+10);
%  tTarg = (round((0.8*foc)/c0/dT):round((1.6*foc)/c0/dT)-1)*dT;
 P = [];
%  tt =554:574;
tt = 553:573;
 Pori = apaz_shifted(:,:,tt,:);
for ii = 1:size(Pori,4)
     P(:,:,:,ii) = permute(interp1(tTarg(tt),permute(squeeze(Pori(:,:,:,ii)),[3 2 1 ]),tTarg(tt(1)):dT/100:tTarg(tt(end)),'spline',0),[3 2 1]);
%     apaz_shifted(:,:,ii) = interp1((t+tShift(ii))',squeeze(apaz_sv(:,:,ii))',tTarg','spline');%,0);
    
end
% P = Pori;
%% project pressure waveforms
dPdz = zeros(size(P));

mask = true(size(dPdz,1));
% calculate the projections using Gmri
% angles = 0:180/nAngles:180-180/nAngles; 
angles = 0:180;
kk = -1:2/(size(dPdz,1)-1):1;
kk = [kron(cosd(angles(:)),kk.') kron(sind(angles(:)),kk.')];
G = Gmri(kk,true(size(dPdz,1)),'fov',size(dPdz,1)/2);
% loop through z-positions and calculated centered differences
for i = 2:length(z_sv)-1
    dPdz(:,:,:,i) = (P(:,:,:,i+1)-P(:,:,:,i-1))/((z_sv(i+1)-z_sv(i-1)));

end
for i = 1:length(z_sv)-2
    for j = 1:size(dPdz,3)
    projdic(:,j,i) = (ifft(ifftshift(G * col(squeeze(dPdz(:,:,j,i+1))))));
%     [xS,info] = qpwls_pcg(0*mask(:),G,1,col(I(:,j,i)),0,0,1,25,mask(:));

%     I = ifft((reshape(I,length(angles),nX)),1);
%     I = I.*abs(kk(:,1)+1i*kk(:,2));
%     r_dPdx(:,:,j,i) = iradon(repmat(squeeze(abs(I')),1,length(angles)),angles,'spline',nY);
    end
end
% projdic = projdic(:,:);
Iz = squeeze(sum(dPdz,2)); % project through the field

%% caculate displacement waveforms in x-z plane
dz = Zd*1/n0*dn0dp*Iz;

dPdx = zeros(size(P));
% loop through x-positions and calculated centered differences
for i = 2:size(P,1)-1
    dPdx(i,:,:,:) = (P(i+1,:,:,:)-P(i-1,:,:,:))/(2*dX);
end
Ix = sum(dPdx,2);
dx = squeeze(Zd*1/n0*dn0dp*sum(dPdx,2));
%% generate 2D histogram at each x-z location
locz = zeros(size(dz));
for i = 2:size(P,4)-1
    locz(:,:,i) = dz(:,:,i)/(z_sv(i+1)-z_sv(i-1));
end
locx = dx/dX;

%%
xxx = round(max(abs(locz(:))));
yyy = round(max(abs(locx(:))));
yyy = xxx;
hisori = zeros(2*yyy+1,2*xxx+1,size(locx,1),size(locx,3));
for kk = 1:size(locx,2)
    for i = 1:size(locx,1)
        for j = 1:size(locx,3)
            x = round(yyy+1+locx(i,kk,j)); z = round(xxx+1+locz(i,kk,j));
            hisori(x,z,i,j)=hisori(x,z,i,j)+1;
        end
    end
end
hisori = hisori;
hissize = 10;
his = zeros(hissize,hissize,size(hisori,3), size(hisori,4)-2);
for i = 1:size(hisori,3)
    for j = 1:size(hisori,4)-2
           tmp = imresize(hisori(:,:,i,j+1)/size(locz,2),[hissize,hissize]);
%           tmp = resample(hisori(:,:,i,j),hissize,size(hisori,1));
%           tmp = resample(tmp',hissize,size(hisori,1));
         tmp = tmp';
        tmp(tmp<0) = 0;
        his(:,:,i,j) = tmp/sum(tmp(:));
    end
end
        
%%
% ghis = zeros(size(his,3)*hissize,size(his,4)*hissize);
% for i = 1:size(locx,1)
%     for j = 1:size(locx,3)
%         ghis(hissize*(i-1)+1:hissize*i,hissize*(j-1)+1:hissize*j) = his(:,:,i,j);
%     end 
% end

 hisdic = his(:,:,:);
 hisdic = premute(hisdic,[3 1 2]);
 hisdic = hisdic(:,:)';
 projdic = projdic(:,:);
% hisdic = permute(hisdic,[3 1 2]);
