% try to use iterations to decrease influences of time delays. It's still being modified...
% projnew(:,:,:,1) = projoritim;
% projprenew(:,:,:,1) = projpretim;
% plot(squeeze(projnew(2,:,71,1))','-');
% hold on;
% pause(0.005);
for i=41:45
    [recon,projnew(:,:,:,i),projprenew,dis(:,:,i),corr(:,:,i),acorr(:,:,:,i)] = iterate_tim(squeeze(projnew(:,:,:,i-1)),squeeze(projprenew));
    plot(squeeze(projnew(2,:,71,i))','-');
    hold on;
    pause(0.005);
end
