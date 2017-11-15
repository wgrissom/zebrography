%example
%hisk(nx*ny*nposition): blurring patterns from captured picetures for slice k.
%% dictonary matching and looking for corresponding projected waveforms.
% clear proj;
% clear drk;
for i = 1:size(his,3)
    [num,proj(:,:,i)] = matchfun(his(:,:,i),hisdic,projdic);
    
end
drk = recon_eachpoint(squeeze(proj(:,153,:)));
