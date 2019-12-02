function proj = forward_model_proj(apaz_sv,dY,dZ,nX,nY,z_sv)

%| Caculate projected pressure by summing up the simulated pressure along the
%| line-of-sight direction.
%|
%| Inputs: 
%| path Saved path of simulated pressure data.
%|
%| Outputs:
%| proj Projected pressure waveforms.

   % load(path,'apaz_sv','dY','dZ','nX','nY','z_sv'); 
    %load simulated pressure by a modified angular specturm method  
    apaz_sv = apaz_sv(2:end-1,:,:);
    nX = nX - 2;
    nY = nY - 2;
    sig = squeeze(max(max(apaz_sv,[],2),[],1));
    nn = find(sig == max(sig));
    % find the slice where the focus is by looking for max peak postive
    % amplitude
    z_svnew = [flip(z_sv(nn):-dZ:z_sv(1)),z_sv(nn)+dZ:dZ:z_sv(end)];
    %make propagation step size uniform
    apaznew = permute(interp1(z_sv,permute(apaz_sv,[3 2 1]),z_svnew,'linear',0),[3 2 1]);
    clear apaz_sv
    [fm,fn] = find(squeeze(squeeze(max(apaznew,[],2))) == max(apaznew(:)));
    sig_focus = squeeze(apaznew(fm,:,fn));
    [~,locs] = findpeaks(sig_focus,'MinPeakHeight',0.95*max(sig_focus(:)));
    nZ = length(z_svnew);
    %find the length of one cycle
    nT = round(mean(diff(locs)));
    nmaxT = find(sig_focus == max(sig_focus(:)));
    apaznew = apaznew(:,round(nmaxT-nT):round(nmaxT+nT),:);
    %interpolation to 3D pressure and calculate projections(in the phase)
    cen = round([nX,nY]/2);
    [y1,x1] = meshgrid(1:nX,1:nX);
    d = (sqrt((x1-cen(1)).^2+(y1-cen(2)).^2));
    r =(-floor(nX/2):floor(nX/2));
    d(d>floor(nX/2)) = -inf;
    delta_nY = nY;
    nproj = nY/delta_nY; %nproj = 1 actually in this study. 
    proj = zeros(nT,nX,nproj,nZ);
    for jj = 1:nZ
        tmp = apaznew(:,round(size(apaznew,2)/2)-round(nT/2):round(size(apaznew,2)/2)+round(nT/2)-1,jj);
        tmp_p = zeros(nT,nX,nY);
        for kk = 1:nT
            tmp_p(kk,:,:,:) = squeeze(interp1(r,squeeze(tmp(:,kk,:)),d,'pchip',0));
        end
        for kk = 1:nproj
            pa = tmp_p(:,:,(kk-1)*delta_nY+1:kk*delta_nY,:);
            proj(:,:,kk,jj) =squeeze(sum(pa(:,:,:)*dY,3));
        end
    end
    proj = single(proj);
