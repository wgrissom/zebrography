function [Pmid,P] = recon_Gmri(apaz_shifted,z_sv,k)

%%
Pori = apaz_shifted;
% [x,y] = meshgrid(1:141,1:141);
% [x2,y2] = meshgrid(1:0.25:141,1:0.25:141);
% for i = 1:size(Pori,3)
%     for j = 1:size(Pori,4)
%         Vq(:,:,i,j) = interp2(x,y,squeeze(Pori(:,:,i,j)),x2,y2);
%     end
% end
zInt = 0:1/5/(size(Pori,4)-1):1;
zInt = zInt(round(length(zInt)/2)-2:round(length(zInt)/2)+2);
%tmp = squeeze(Pori(71,71,:,:));
%Pupsamp = interp1(0:1/(size(tmp,2)-1):1,permute(tmp,[2 1]),zInt,'spline');
Pupsamp = interp1(0:1/(size(Pori,4)-1):1,permute(Pori,[4 1 2 3]),zInt,'spline');
Pupsamp = permute(Pupsamp,[2 3 4 1]);
P = Pupsamp;

dPdx = zeros(size(P));
nt = size(P,3);
nY = size(P,2);
nX = size(P,1);
% loop through z-positions and calculated centered differences
for i = 2:size(P,4)-1
    dPdx(:,:,:,i) = (P(:,:,:,i+1)-P(:,:,:,i-1))/((z_sv(i+1)-z_sv(i-1)))*5;
end
% dPdx = dPdx(:,:,:,2:4); % select middle three slices of differences
% I = squeeze(sum(dPdx,2)); % project through the field
nAngles = ceil(pi/2*nY);

mask = true(size(dPdx,1));
% calculate the projections using Gmri
% angles = 0:180/nAngles:180-180/nAngles; 
angles = 0:180/nAngles:180 - 180/nAngles;
kk = -1:2/(size(dPdx,1)-1):1;
kk = [kron(cosd(angles(:)),kk.') kron(sind(angles(:)),kk.')];
G = Gmri(kk,true(size(dPdx,1)),'fov',size(dPdx,1));
% angles = 0:180; 
% kk = -1:2/(size(I,1)-1):1;
% kk = [kron(cosd(angles(:)),kk.') kron(sind(angles(:)),kk.')];
% G = Gmri(kk,true(size(I,1)),'fov',size(I,1)/2);
for i = 2:size(P,4)-1
    for j = 1:nt
    I = G * col(squeeze(dPdx(:,:,j,i)));
%     [xS,info] = qpwls_pcg(0*mask(:),G,1,col(I(:,j,i)),0,0,1,25,mask(:));
    I = I.*abs(kk(:,1)+1i*kk(:,2));
%     I = ifft((reshape(I,length(angles),nX)),1);
%     I = I.*abs(kk(:,1)+1i*kk(:,2));
%     r_dPdx(:,:,j,i) = iradon(repmat(squeeze(abs(I')),1,length(angles)),angles,'spline',nY);
    r_dPdx(:,:,j,i) = real(reshape(G'*I,nX,nY)/(nX*nY));
    end
end
% to get projections:
% 1) multiply your image into G. This gives you a vector of all the
% projections' Fourier transforms
% 2) take inverse FT of each projection, with fftshifts

% to get back reconstructed images:
% 1) Take forward FT of each projection
% 2) Multiply by |k| filter for density compensation
% 3) Multiply stacked projections by G'
% 4) reshape


%% muti-timepoints
% r_dPdx =zeros(size(P));
% for j = 1:nt
%     for i = 2:size(P,4)-1
%         r_dPdx(:,:,j,i) = iradon(repmat(squeeze(I(:,j,i)),1,nAngles),180/nAngles,'spline',nY);
%         % r_dPdx(:,:,j,i) = iradon(repmat(squeeze(I(:,j,i)),1,nAngles),0:180/nAngles:180 - 180/nAngles,nX);
%         %outputsize = nX;
%     end
% end
% Pmid = zeros(size(P));
% for i = 3:size(P,4)-2
%     Pmid(:,:,:,i) = -squeeze((1/k)^2 * (r_dPdx(:,:,:,i+1) - r_dPdx(:,:,:,i-1))/((z_sv(i+1)-z_sv(i-1))));
% end
Pmid = -squeeze((1/k)^2 * (r_dPdx(:,:,:,4) - r_dPdx(:,:,:,2))/(z_sv(4)-z_sv(2)))*5;
P = P(:,:,:,3);
% figure; plot( squeeze(P(round(nX/2),round(nY/2),:,3)));
% hold on;
% plot( squeeze(Pmid(round(nX/2),round(nY/2),:)),'r');
% xlabel('tTarg');
% ylabel('Pressure');
% legend('real pressure','reconstructed pressure');
%
%
%  %% At the time instant when pressure wave reaches the focus.
% sc = squeeze(P(round(nX/2),round(nY/2),:,3));
% ntf = find(sc == max(sc));
% figure; imshow(P(:,:,ntf,3),[]); colorbar;
% title('Real pressure');
% figure; imshow(Pmid(:,:,ntf),[]); colorbar;
% title('Reconstructed pressure');


%%c