% 
% [hisdic,projdic] = generate_dict(apaz_shifted,z_sv,dX,dT,tTarg);
hisColl = hisdic(:,:,:);
hisColl = permute(hisColl,[3 1 2]);
hisColl = hisColl(:,:)';
projColl = permute(projdic, [2 1 3]);
projColl = projColl(:,:);
%%
histest = hisdic(:,:,:,424:564);
quarter = round(1000/4);
nt = 1000;
Phase = round(quarter*round(rand(1,size(projColl,2))));% - 20 + 40 *rand(1,size(projColl,2)));
zerox = rand(1,round(size(projColl,2)*9/10));
zerox = round(zerox*size(projColl,2));
zerox = zerox +1;
Phase(zerox) = 0;
for i = 1:size(projColl,2)
    tmp = projColl(:,i);
    tmp = [tmp;tmp];
    tmp = real(ifft(fft(tmp).*exp(-1i*Phase(i)*(0:(length(tmp)-1))'*2*pi/length(tmp))));
    projtimdic(:,i) = tmp(nt+1:end);
end
for i = 1:141
%     display(num2str(i));
    for j = 1:180
        [num,projoritim(j,:,i)] = matchfun(histest(:,:,j,i),hisColl,projtimdic);
%         display([num2str(i),',',num2str(j)]);
    end
    
end

for i = 1:nt
    recontim(:,:,i)= recon_eachpoint(squeeze(projoritim(:,i,:)));
    % figure; im(drk);
    % figure; im(P(:,:,153,12));
end

P = recontim;
angles = 0:179;
kk = -1:2/(size(P,1)-1):1;
kk = [kron(cosd(angles(:)),kk.') kron(sind(angles(:)),kk.')];
G = Gmri(kk,true(size(P,1)),'fov',size(P,1)/2);
for j = 1:size(P,3)
    I = reshape(G * (squeeze(P(:,:,j))),[size(P,1),size(angles)]);
    projpretim(:,j,:) = squeeze(real(ifftshift(ifft(ifftshift(I,1),[],1))))';   
end
save('utimedelay','-v7.3')
