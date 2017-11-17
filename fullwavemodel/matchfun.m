function [num,proj] = matchfun(his,hisdic,projdic)
% hisdic: historgrams in dictionaries nx*ny*number of historgrams
% projdic: projected waveforms
hisColl = permute(hisdic,[3 1 2]);
hisColl = hisColl(:,:)';
num = 0;
if length(find(his~=0)) == 1
    proj = zeros(size(projdic,1),1);
else
    hasPressure = sqrt(sum((hisColl - repmat(col(his),[1 size(hisColl,2)])).^2,1));
    num = find(hasPressure == min(hasPressure));
    proj = squeeze(projdic(:,num(1)));
end