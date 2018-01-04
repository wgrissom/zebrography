% [hisdic,projdic] = generate_dict(apaz_shifted,z_sv,dX,dT,tTarg);
hisColl = hisdic(:,:,:);
hisColl = permute(hisColl,[3 1 2]);
hisColl = hisColl(:,:)';
projColl = permute(projdic, [2 1 3]);
projColl = projColl(:,:);
histest = hisdic(:,:,:,424:564);
for i = 1:141 
    display(num2str(i));
    for j = 1:180
        [num,projori(j,:,i)] = matchfun(histest(:,:,j,i),hisColl,projColl);
         display([num2str(i),',',num2str(j)]);
    end

end

for i = 1:size(projdic,2)
    reconori(:,:,i)= recon_eachpoint(squeeze(projori(:,i,:)));                                                                 
% figure; im(drk);
% figure; im(P(:,:,153,12));
end

P = reconori;
angles = 0:179;
kk = -1:2/(size(P,1)-1):1;
kk = [kron(cosd(angles(:)),kk.') kron(sind(angles(:)),kk.')];
G = Gmri(kk,true(size(P,1)),'fov',size(P,1)/2);
for j = 1:size(P,3)
    I = reshape(G * (squeeze(P(:,:,j))),[size(P,1),size(angles)]);
    projpreori(:,j,:) = squeeze(real(ifftshift(ifft(ifftshift(I,1),[],1))))';   
end

