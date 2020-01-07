function [res,thetahat,coeff,A] = l0polyfit(theta, img, order, mask, ...
                                            nrwt, epsl1, unwrapswitch)

%
% [res,thetahat,coeff,A] = l0polyfit(theta, img, order, mask, nrwt, epsl1, unwrapswitch)
%
% Performs reweighted-l1 polynomial fitting
%
% Inputs:  
%
% theta - phase image
% img - magnitude image
% order - polynomial order
% mask - error mask, also used to calculate sigma
% nrwt - number of times to reweight the l1 problem
%        default is 25
% epsl1 - reweighted l1 damping parameter
%         default is 0.01
% unwrapswitch - perform phase unwrapping prior to 
%                estimation
%
% Outputs:
% 
% res - residual difference between theta and estimated 
%       polynomial
% thetahat - estimated polynomial
% coeff - estimated polynomial coefficients
% A - polynomial matrix

if nargin < 4 
  mask = ones(size(img));
end
if nargin < 5
  nrwt = 25;
end
if nargin < 6
  epsl1 = 0.01;
end
if nargin < 7
  unwrapswitch = 0;
end

if unwrapswitch
  theta = run_unwrap2(img.*exp(1i*theta));
end

% construct polynomial matrix
[sx,sy] = size(theta);
[yc,xc] = meshgrid(linspace(-1/2,1/2,sy), linspace(-1/2,1/2,sx));
yc = yc(:);
xc = xc(:);
A = [];
if sy > 1 % if not a one-d fit
  for yp = 0:order
    for xp = 0:(order-yp)
      A = [A (xc.^xp).*(yc.^yp)];
    end
  end
else
  for xp = 0:order
    A = [A (xc.^xp)];
  end
end  

npoly = size(A,2);
disp(sprintf('Fitting %d polynomial basis functions',npoly));

% combine image magnitude and mask weighting
w = abs(img).*mask;

% apply mask and image magnitude weighting
WA = A.*repmat(w(:),[1,size(A,2)]);
Wtheta = theta(:).*w(:);

for ii = 1:nrwt+1
  
  % using cvx to solve l1 problem is slow
  
  %cvx_begin
  % variable coeff(npoly);
  % minimize( norm( WA*coeff - Wtheta(:),1) );
  %cvx_end

  % get initial coeffs, weights
  coeff = inv(WA'*WA)*WA'*Wtheta;
  Wthetahat = WA*coeff;
  rl2 = Wtheta - Wthetahat;
  
  % get 'step size'
  epsl2_init = std(rl2);
  epsl2 = epsl2_init;
  
  % start iterations
  while epsl2 > (10^-15 * epsl2_init)
    rl2_old = 2 * rl2; % just a bs number to get the loop going
    while (norm(rl2_old)-norm(rl2)) > sqrt(epsl2)/100*norm(rl2_old);
      rl2_old = rl2;
      % calculate weights from previous residual.
      % one sqrt comes from the norm we approx, the other comes
      % from the fact that we multiply this twice, but only want
      % it weighted once
      wtsl2 = sqrt(1./(epsl2 + abs(rl2).^2).^(1/2));
      WWA = WA.*repmat(wtsl2,[1,npoly]);
      % calculate coeffs
      coeff = inv(WWA'*WWA)*WWA'*(Wtheta.*(wtsl2));
      % calculate new residual
      Wthetahat = WA*coeff;
      rl2 = Wtheta - Wthetahat;
    end
    % decrease epsilon
    epsl2 = epsl2/10;
  end

  % calculate weighted residual
  thetahat = A*coeff;
  rl1 = w(:).*abs(theta(:) - thetahat);
  
  % calculate weights for next iteration
  wtsl1 = w(:)./(rl1 + epsl1);
  %im(reshape(wtsl1,[64 64]));drawnow
  
  % apply mask, image magnitude, and reweighting
  WA = A.*repmat(wtsl1,[1,npoly]);
  Wtheta = theta(:).*wtsl1;
  
end

thetahat = reshape(thetahat,[sx sy]);

res = theta - thetahat;

