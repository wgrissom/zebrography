function [dxreal,dzreal] = forward_model_dxdz(apaz_sv,dX,dY,dZ,nX,nY,z_sv,Zd)

%| Caculate projected pressure by summing up the simulated pressure along the
%| line-of-sight direction.
%|
%| Inputs:
%| path   Saved path of simulated pressure data.
%| Zd     The distance from the middle of the FUS beam to the iPad screen. 
%|
%| Outputs:
%| apaznew
%| dxreal,dzreal Projected pressure waveforms.
  %  load(path,'apaz_sv','dX','dY','dZ','nX','nY','z_sv','dT','nT')
    apaz_sv = apaz_sv(2:end-1,:,:);
    nX = nX - 2;
    nY = nY - 2;
    sig = squeeze(max(max(apaz_sv,[],2),[],1));
    nn = find(sig == max(sig));
    display(nn); % find the slice where the focus is.(Peak postive
    z_svnew = [flip(z_sv(nn):-dZ:z_sv(1)),z_sv(nn)+dZ:dZ:z_sv(end)];%make propagation step size uniform
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
    delta_nY = 5;
    nproj = nY/delta_nY;
    delta_y = Zd/nproj;
    n = 4;
    dx = zeros(nT,nX,nproj,nZ);
    dz = zeros(nT,nX,nproj,nZ);
    for jj = 1+n/2:nZ-n/2
        tmp = apaznew(:,round(size(apaznew,2)/2)-round(nT/2):round(size(apaznew,2)/2)+round(nT/2)-1,jj-n/2:jj+n/2);
        tmp_p = zeros(nT,nX,nY,n+1);
        for kk = 1:nT
            tmp_p(kk,:,:,:) = squeeze(interp1(r,squeeze(tmp(:,kk,:)),d,'pchip',0));
        end
        for kk = 1:nproj
            pa = tmp_p(:,:,(kk-1)*delta_nY+1:kk*delta_nY,:);
            [dxp,dzp] = accum_d(pa,dX,dY,dZ,n);
            dx(:,:,kk,jj) = dxp;
            dz(:,:,kk,jj) = dzp;
        end
    end
    apaznew = single(apaznew);
    [z0,x0] = meshgrid(z_svnew,(-round(nX/2):round(nX/2))*dX);
    x0 = x0(2:end-1,:);
    z0 = z0(2:end-1,:);
    nXp = size(z0,1);
    nZp = size(z0,2);
    dxx = zeros(nT,nXp,nZp,nproj);
    dzz = zeros(nT,nXp,nZp,nproj);
    dxreal = zeros(nT,nXp,nZp);
    dzreal = zeros(nT,nXp,nZp);
    for kk = 1:nT
        cenZp = z0;
        cenXp = x0;
        for nn = 1:nproj
            dxp = interp2(z0,x0,squeeze(dx(kk,:,nn,:)),cenZp,cenXp,'cubic',0);
            dzp = interp2(z0,x0,squeeze(dz(kk,:,nn,:)),cenZp,cenXp,'cubic',0);
            dxx(kk,:,:,nn) = dxp;
            dzz(kk,:,:,nn) = dzp;
            dxtmp = delta_y*nn*(sum(dxx(kk,:,:,1:nn),4));
            dztmp = delta_y*nn*(sum(dzz(kk,:,:,1:nn),4));
            cenXp = x0 + squeeze(dxtmp);
            cenZp = z0 + squeeze(dztmp);
        end
        dxreal(kk,:,:) =(dxtmp);
        dzreal(kk,:,:) =(dztmp);
    end
    dxreal = single(dxreal);
    dzreal = single(dzreal);

%     save([Basepath,'/params_',num2str(P0(pp)),'_',num2str(fnum(bb)),'.mat'],'dX','dY','dZ','nX','nZ','nY','p0','-v7.3');
%     save([Basepath,'/pressure_',num2str(P0(pp)),'_',num2str(fnum(bb)),'.mat'],'apaznew','-v7.3');
%     save([Basepath,'/displacement_',num2str(P0(pp)),'_',num2str(fnum(bb)),'.mat'],'dxreal','dzreal','-v7.3');
