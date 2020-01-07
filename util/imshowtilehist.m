function histmp = imshowtilehist(his)
%| function imshowtilehist
%|
%| Input:
%| his    [nblock,nblock,nx,nz] sets of segmented histograms
%| 
%| Output:
%| tiled histogram  [nblock*nx,nblock*nz]
histmp =  permute(his,[1 3 2 4]);
histmp = histmp(:,:,:);
histmp = permute(histmp,[3 1 2]);
histmp = histmp(:,:)';
