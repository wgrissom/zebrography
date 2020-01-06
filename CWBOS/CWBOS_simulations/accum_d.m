function [dx,dz] = accum_d(pa,dX,dY,dZ,n)
%| Calculate physical displacements of pixels in the background pattern.
%|
%| Inputs:
%|        pa, [Nt,Nx,Ny,Nz] simulated pressure
%|        dX, dY  grid spacing
%|        dZ  wave propagation's step size
%|        n   n-th order finite difference
%| 
%| Outputs: 
%|        dx [Nt,Nx] physical displacements(meter) in the x-dimension.
%|        dz [Nt,Nx] physical displacements(meter) in the z-dimension.

%refractive angles
tmp_proj =squeeze(sum(pa*dY,3));
nX = size(pa,2);
A = zeros(n+1,n+1);
for ii = 1:n+1
    A(ii,:) = (-n/2:n/2).^(ii-1);
end
c = A\[0;1;zeros(n-1,1)];

n0 = 1.3325; % base index of refraction

dn0dp = 14.83/101325*1e-6;% dindex/dPascal, from Waxler paper, read from Figure 2 at 24.8 deg CSS
for kk = 1:n+1
    tmp_proj(:,:,kk) = circshift(tmp_proj(:,:,kk),[10*(kk-n/2-1) 0]);
end

dx = [];
dz = [];
nT = size(tmp_proj,1);
dx = zeros(nT,nX);
dz = zeros(nT,nX);
for ii = 1+n/2:nX-n/2
    dx(:,ii) = 1/n0*dn0dp*squeeze(tmp_proj(:,ii-n/2:ii+n/2,n/2+1))*c/dX;
    dz(:,ii) = 1/n0*dn0dp*squeeze(tmp_proj(:,ii,:))*c/dZ;
end
