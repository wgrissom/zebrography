apaz = apaz_sv(:,:,2:end);
z_sv = cumsum(zvec);
z_sv = z_sv(abs(z_sv-foc) < 0.25*foc);
num = find(abs(z_sv-foc) == min(abs(z_sv-foc)));
tTarg = (round((0.5*foc)/c0/dT):round((2*foc)/c0/dT)-1)*dT;

%% shift in z
apaz = reshape(apaz,[587,1,2881,116]);
zPos = z_sv;
tShift = zPos/c0; % time shifts (position/(m/s) = s)
apaz(isnan(apaz)) = 0;

%%
apaz_shifted = [];
for ii = 1:size(apaz,4)
     apaz_shifted(:,:,:,ii) = permute(interp1((t+tShift(ii))',permute((apaz(:,:,:,ii)),[3 2 1 ]),tTarg','spline',0),[3 2 1]);
%     apaz_shifted(:,:,ii) = interp1((t+tShift(ii))',squeeze(apaz_sv(:,:,ii))',tTarg','spline');%,0);
    
end

%%
n0 = 1.3325; % base index of refraction
% m, distance from center of US focus to background
Zd = 0.05;
dn0dp = (1.3475-1.3325)/((1100-1)*100000); % dindex/dPascal, from Waxler paper, read from Figure 2 at 24.8 deg CSS
%%
P = [];
tt = 1911:1953;
Pori = apaz_shifted(:,:,tt,:);
% P = Pori;
for ii = 1:size(Pori,4)
    P(:,:,:,ii) = permute(interp1(tTarg(tt),permute((Pori(:,:,:,ii)),[3 2 1]),tTarg(tt(1)):dT/50:tTarg(tt(end)),'PCHIP',0),[3 2 1]);%   apaz_shifted(:,:,ii) = interp1((t+tShift(ii))',squeeze(apaz_sv(:,:,ii))',tTarg','spline');%,0);
end
[x,y] = findpeaks(squeeze(P(round(size(P,1)/2),1,:,num)));
P = P(:,:,y(1):y(end),:);

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