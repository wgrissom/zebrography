%example 
%hisk(nx*ny*nposition): blurring patterns from captured picetures for slice k.
%% dictonary matching and looking for corresponding projected waveforms.
for i = 1:size(his,3)
    [num,proj(:,:,i)] = matchfun(his(:,:,i),hisdic,projdic);
    for j = 1:size(proj,2)
     drk(:,:,j,i) = recon_eachpoint(proj(:,j,i));
    end
end