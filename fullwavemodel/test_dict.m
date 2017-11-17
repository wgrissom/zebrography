%example
%hisk(nx*ny*nposition): blurring patterns from captured picetures for slice k.
% %% dictonary matching and looking for corresponding projected waveforms.
clear proj;
clear drk;
for i = 1:141
    for j = 1:180
    [num,proj(j,:,i)] = matchfun(histest(:,:,j,i),squeeze(hisdic(:,:,j,:)),squeeze((projdic(j,:,:))));
    end
end
drk = recon_eachpoint(squeeze(proj(:,153,:)));
figure; im(drk);
figure; im(P(:,:,153,12));


% %%
% for i = 1:size(his,3)
%     hisColl = permute(hisdic,[3 1 2]);
%     hisColl = hisColl(:,:)';
%     hasPressure = sqrt(sum((hisColl - repmat(col(his(:,:,i)),[1 size(hisColl,2)])).^2,1));
%     num = find(hasPressure == min(hasPressure));
%     if length(num)>1
%         proj1(:,:,i) = -1e-7;
%     else
%         proj1(:,:,i) = squeeze(projdic(:,:,num));
%     end
% end
% drk = recon_eachpoint(squeeze(proj1(:,153,:)));
% figure; im(drk);
% figure; im(squeeze(P(:,:,153,12)));
