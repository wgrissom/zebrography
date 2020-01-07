function [img,epsilonPix, d] = calcImgMC(p,alpha, background)

% Overal loop (more general than below since any background image can be used):
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For each spatial location:
%   Calculate displacements
%   For each line position/output image:
%       Add displaced pixels into image
% 3rd dim: line position
epsilonPix = -alpha*squeeze(diff(diff(p,1,1),1,2));
[Nx,Ny,Nt] = size(epsilonPix);
d = epsilonPix;
tmp = background';
bg = background';
for ii = 1:Nx % Loop over line positions
    d(ii,:,:) = round(d(ii,:,:) + ii);
    for jj = 1:Ny
        for tt = 1:Nt
    % calculate sinusoid
    %s = sin(2*pi*f*t);% - k*x(ii)); % k doesnt matter here since we integrate over so many cycles
    %c = cos(omega*t);% - k*x(ii));
     
    %for kk = 1:Nx % Loop over vertical positions
        
        % calculate displacements
   % displacements - was a dt here - why?
        %epsilon = epsilon + 1/n0*c*Zd*(p(ii+1,kk)-p(ii-1,kk))/(2*dx)*dn0dp;
        %epsilonPix = Zd*epsilon/dx; % displacement in pixels
        
        % average displacements

        if d(ii,jj,tt) > Nx
            d(ii,jj,tt) = Nx;
        end
        if d(ii,jj,tt) < 1
            d(ii,jj,tt) = 1;
        end
        tmp(d(ii,jj,tt),jj) = tmp(d(ii,jj,tt),jj) + bg(ii,jj);
        end
    end
end
      
img = tmp';
