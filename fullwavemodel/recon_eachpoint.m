function drk = recon_eachpoint(proj)
%reconstruct pressure map of slice k using GMRI
%proj(t,x) timepoints* # rows
angles = 0:179;
nx = (size(proj,2));
% nx = size(proj,2)/length(angles);
kk = -1:2/(nx-1):1;
kk = [kron(cosd(angles(:)),kk.') kron(sind(angles(:)),kk.')];
G = Gmri(kk,true(nx,nx),'fov',[nx/2,nx/2]);
I = col(fftshift(fft(fftshift(proj'),[],1),1));
I = I.*abs(kk(:,1)+1i*kk(:,2));
drk = real(reshape(G'*I,[nx,nx])/(nx*nx));
                                                                                                                                     