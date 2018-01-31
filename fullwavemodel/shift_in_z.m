function apaz_shifted = shift_in_z(apaz_sv,z_sv,tTarg,t,c0)

%load adaptive_results_orig_zQuarterSteps_quarterPressure
% contains: apaz_sv: Nx x Nz x Nt

% for each z-location, shift the time dimension. The curves should 
% move to the right (forward in time) as z increases
% for each z-location, shift the time dimension. The curves should 
% move to the right (forward in time) as z increases

%  dX=lambda/4; dY=lambda/4;% grid spacing

zPos = z_sv;
tShift = zPos/c0; % time shifts (position/(m/s) = s)
apaz_sv(isnan(apaz_sv)) = 0;

%%
apaz_shifted = [];
for ii = 1:size(apaz_sv,4)
     apaz_shifted(:,:,:,ii) = permute(interp1((t+tShift(ii))',permute(squeeze(apaz_sv(:,:,:,ii)),[3 2 1 ]),tTarg','pchip',0),[3 2 1]);
%     apaz_shifted(:,:,ii) = interp1((t+tShift(ii))',squeeze(apaz_sv(:,:,ii))',tTarg','spline');%,0);
    
end
