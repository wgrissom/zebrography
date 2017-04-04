disp 'Loading images'
% load the displayed image
img_disp = imread('hex1.png');
% get the fraction of dark pixels, for setting a threshold below
pinFrac = sum(img_disp(:) == 0)/numel(img_disp);


% load the measured image
img_meas = imread('IMG_0064.CR2');
img_meas = mean(img_meas(50:3800, 150:end, :),3); % crop black edges, combine color channels

% crop it to middle
[sx,sy] = size(img_meas);
img_meas = img_meas(floor(sx/2)-round(sx/4):floor(sx/2)+round(sx/4),...
   floor(sy/2)-round(sy/4):floor(sy/2)+round(sy/4));

skip = 100; % # voxels to skip for low res ops like poly fit

disp 'Doing polynomial fit'
% remove slow intensity modulations
order = 3; % poly order
[~,mods,c] = l0polyfit(img_meas(1:skip:end,1:skip:end),...
    ones(size(img_meas(1:skip:end,1:skip:end))),order,...
    ones(size(img_meas(1:skip:end,1:skip:end))));
[sx,sy] = size(img_meas)
[yc,xc] = meshgrid(linspace(-1/2,1/2,sy), linspace(-1/2,1/2,sx));
% yc = yc(:);
% xc = xc(:);
img_meas_remmod = img_meas; 
ind = 1;
    for yp = 0:order
      for xp = 0:(order-yp)
        img_meas_remmod = img_meas_remmod - c(ind)*(xc.^xp).*(yc.^yp);
        ind = ind + 1;
      end
    end
 
 
% get a threshold
tmp = sort(img_meas_remmod(:));
thresh = tmp(round(pinFrac*length(tmp)));
%pinMaskInit = img_meas_remmod < thresh;

disp 'Finding rotation angle'
% get the angle of the pattern by rotating the mask itself, until # of zeros is maximized. 
angles = -2.5:0.1:2.5;
numZeros = zeros(size(angles));
for ii = 1:length(angles)
    img = imrotate(img_meas_remmod,angles(ii));
    pinMaskTmp = img < thresh;
    rowCollapsed = sum(pinMaskTmp,1); % collapse rows
    colCollapsed = sum(pinMaskTmp,2); % collapse columns
    numZeros(ii) = sum(rowCollapsed == 0) + sum(colCollapsed == 0);
end
gridAngle = angles(numZeros == max(numZeros));
imgAligned = imrotate(img_meas_remmod,gridAngle,'bilinear');
pinMaskAligned = imgAligned < thresh;
imgMaskedAligned = pinMaskAligned.*imgAligned;

disp 'Calculating grid spacing'
% get grid spacing in each dimension
pinMaskRowCollapsed = sum(pinMaskAligned,1); % collapse rows
pinMaskColCollapsed = sum(pinMaskAligned,2); % collapse columns
spacing = 5:0.05:length(pinMaskRowCollapsed)/10;
for ii = 1:length(spacing)
    pinMaskCollapsedRowSynth = false(size(pinMaskRowCollapsed));
    pinMaskCollapsedRowSynth(round(1:spacing(ii):length(pinMaskRowCollapsed))) = true;
    pinMaskCollapsedColSynth = false(size(pinMaskColCollapsed));
    pinMaskCollapsedColSynth(round(1:spacing(ii):length(pinMaskColCollapsed))) = true;
    spacingCorr(ii) = abs(ft(pinMaskRowCollapsed))*abs(ft(pinMaskCollapsedRowSynth))'./norm(double(pinMaskCollapsedRowSynth)) + ...
        abs(ft(pinMaskColCollapsed))'*abs(ft(pinMaskCollapsedColSynth))./norm(double(pinMaskCollapsedColSynth));
end
dCol = spacing(spacingCorr == max(spacingCorr));

% another idea: generate pin pattern with detected spacing, and find
% shift of that pattern wrt mask by absorbing FT phase or by collapsing and
% shifting. Optional: loop through the pin centers and find centroids? 
pinMaskSynth = false(size(pinMaskAligned));
pinMaskSynth(round(1:2*dCol:size(pinMaskSynth,1)),round(1:2*dCol:size(pinMaskSynth,2))) = true;
pinMaskSynth(round(dCol + (1:2*dCol:size(pinMaskSynth,1)-dCol)),round(dCol + (1:2*dCol:size(pinMaskSynth,2)-dCol))) = true;
nPins = sum(pinMaskSynth(:));
pinWidth = round(sqrt(pinFrac*numel(pinMaskSynth)/nPins));
pinMaskSynthPSF = conv2(double(pinMaskSynth),ones(pinWidth),'same');
pinMaskSynthRowCollapsed = sum(pinMaskSynthPSF,1);
pinMaskSynthColCollapsed = sum(pinMaskSynthPSF,2);

%tmp = ift2(abs(ft2(pinMaskSynth)).*exp(1i*angle(ft2(pinMaskAligned))));
shift = 0:floor(dCol);
for ii = 1:length(shift)
    shiftCorrCol(ii) = pinMaskRowCollapsed*...
        circshift(pinMaskSynthRowCollapsed,shift(ii),2)';
    shiftCorrRow(ii) = pinMaskColCollapsed'*...
        circshift(pinMaskSynthColCollapsed,shift(ii),1);
end
pinMaskSynthAligned = circshift(pinMaskSynth,...
    [shift(shiftCorrRow == max(shiftCorrRow)),...
    shift(shiftCorrCol == max(shiftCorrCol))]);
pinMaskSynthAligned(1:shift(shiftCorrRow == max(shiftCorrRow)),:) = 0;
pinMaskSynthAligned(:,1:shift(shiftCorrCol == max(shiftCorrCol))) = 0;
pinMaskSynthPSFAligned = circshift(pinMaskSynthPSF,...
    [shift(shiftCorrRow == max(shiftCorrRow)),...
    shift(shiftCorrCol == max(shiftCorrCol))]);
pinMaskSynthPSFAligned(1:shift(shiftCorrRow == max(shiftCorrRow)),:) = 0;
pinMaskSynthPSFAligned(:,1:shift(shiftCorrCol == max(shiftCorrCol))) = 0;

figure;imagesc(pinMaskAligned+pinMaskSynthPSFAligned);axis image;colormap gray
title 'Aligned + Synthesized'
figure;imagesc(pinMaskAligned);axis image;colormap gray
hold on
[x,y] = ndgrid(1:size(pinMaskAligned,1),1:size(pinMaskAligned,2));
plot(y(pinMaskSynthAligned == 1),x(pinMaskSynthAligned == 1),'*g');


% Still todo:
% - convert pin locations to coordinates
% - rotate pin coordinates back to same angle as original image
