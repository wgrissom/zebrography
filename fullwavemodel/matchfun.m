function [num,proj] = matchfun(his,hisdic,projdic)
% hisdic: historgrams in dictionaries nx*ny*number of historgrams
% projdic: projected waveforms
hisColl = permute(hisdic,[3 1 2]);
hisColl = hisColl(:,:)';
hasPressure = sqrt(sum((hisColl - repmat(col(his),[1 size(hisColl,2)])).^2,1));
num = find(hasPressure == min(hasPressure));
proj = squeeze(projdic(:,:,num));