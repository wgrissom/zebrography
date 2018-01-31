
%%% Interpolation to match with experimental image. 
dgrid = 92.6/98*1e-3;
nn = find(abs(z_sv-2*a) == min(abs(z_sv-2*a)));
z_svnew = [flip(z_sv(nn):-dgrid:z_sv(1)),z_sv(nn)+dgrid:dgrid:z_sv(end)];
apaz = permute(interp1(0:dX:(nX*dX-dX),permute(squeeze(apaz_shifted),[1 2 3 4]),[flip(69*dX:-dgrid:0),69*dX+dgrid:dgrid:(nX-1)*dX],'PCHIP',0),[1 2 3 4]);
apaz = permute(interp1(0:dX:(nX*dX-dX),permute(squeeze(apaz),[2 1 3 4]),[flip(69*dX:-dgrid:0),69*dX+dgrid:dgrid:(nX-1)*dX],'PCHIP',0),[2 1 3 4]);
apaz = permute(interp1(z_sv,permute(squeeze(apaz),[4 3 2 1]),z_svnew,'PCHIP',0),[4 3 2 1]);
[img,P] = generate_simu(apaz,z_svnew,dgrid,dT,tTarg,15);