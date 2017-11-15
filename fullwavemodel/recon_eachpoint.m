function drk = recon_eachpoint(proj)
%reconstruct pressure map of slice k using GMRI
%proj(t,x) timepoints* # rows
angles = 0:180;
nx = length(size(proj,1))/length(angles);
% nx = size(proj,2)/length(angles);
kk = -1:2/(nx-1):1;
kk = [kron(cosd(angles(:)),kk.') kron(sind(angles(:)),kk.')];
G = Gmri(kk,true(nx),'fov',nx/2);
I = fftshift(fft(squeeze(proj)));
I = I.*abs(kk(:,1)+1i*kk(:,2));
drk = real(reshape(G'*I,[nx,nx])/(nx*nx));
                                                                                                                                     