%|Process the actual photos acquired in the experiments. 
%|
%|
%%This is a simple simulated non-FUS histograms used for actual photo
%%processing. Size of pins and spacing between pins were determined from the actual
%%photo.
nblock = 54; % width of histograms, even number to be consistent with ipad and camera
npin = 6;% pin: 1 pixel, spacing: 8 pixels
blocks = ones(nblock,nblock);
nhalf = 27;
blocks((nblock-npin+2)/2:(nblock+npin)/2,(nblock-npin+2)/2:(nblock+npin)/2) = 0;

photopath ='/200mv'; % Path where the photos saved.
list1 = dir([photopath,'/*.CR2']);
datapath = [photopath,'/data'];
mkdir(datapath);

repnum = 1; %number of repeated photos.
circnum = 4; % 4x4 photos with equal-interval grid translation in x- and z- dimension.

order = reshape(1:circnum^2,circnum,circnum)'; % The order of photo acquisition. This is z-direction first and then x-direction.
order = order(:);

% x-direction first and then z-direction.
% order = reshape(1:circnum^2,circnum,circnum);
% order = order(:);


%import one pair of nonFUS and FUS photos to determine FOV manually.
bkg =double(rgb2gray(imread([photopath,'/',list1(1).name])));  %Load non-FUS photo, (0,0) grid-translation
blur =double(rgb2gray(imread([photopath,'/',list1(2).name]))); %Load FUS photo, (0,0) grid-translation

gen_FOV = false;

xxmin = 1479; xxmax = 2167; zzmin = 2292; zzmax = 3660; %An example of extracted FOV for 200mV


%-----------------------------------------------------------%
% Use 'roipoly' function to extract a rectangular FOV (pin-blurring region)
% around the focus.
if gen_FOV
    imshow(blur,[]);
    [~,zind,xind] = roipoly;
    zind = round(zind); xind = round(xind);
    xxmax = max(xind(:)); xxmin = min(xind(:));
    zzmax = max(zind(:)); zzmin = min(zind(:));
    close all;
end
%--------------------------------------------------------%

bkg = bkg(xxmin:xxmax,zzmin:zzmax); %Crop the full FOV into extracted FOV. 
blur = blur(xxmin:xxmax,zzmin:zzmax);


%%Preprocessing non-FUS photo to calculate the actual size of histogram in the acquired photo. 
bkgtmp = filter2(fspecial('average',3),bkg/max(bkg(:)));
bkgtmp = (bkgtmp-min(bkgtmp(:)))/(max(bkgtmp(:))-min(bkgtmp(:)));
[N,X] = hist(bkgtmp(:));
edges = filter2(fspecial('average',5),(bwmorph(bkgtmp>median(X),'shrink',inf) == 1))>0.1;  
%% smooth the original photo and use the Matlab function "bwmorph" to extract the edges of histograms.
[m,n] = find(edges == 1);
tbl1 = tabulate(m);
tbl2 = tabulate(n);
thres1 = smoothdata(tbl1(:,2));
thres2 = smoothdata(tbl2(:,2));
thres1 = thres1/max(thres1(:));
thres2 = thres2/max(thres2(:));
[~,zz] = findpeaks(thres2,'MinPeakHeight',0.7);
[~,xx] = findpeaks(thres1,'MinPeakHeight',0.7);  % find the coordinates (x,z) of endpoints.
diffzz = diff(zz(2:end-1));
diffxx = diff(xx(2:end-1));
nblock_p = round((mean(diffzz)+mean(diffxx))/2);
if mod(nblock_p,2)~=0
    nblock_p = nblock_p+1; % This is the true width of simulated histograms. 
end
nxp = size(bkg,1); % the number of histograms in x-dimension
nzp = size(bkg,2); % the number of histograms in z-dimension
[x1,z1] = meshgrid(1:nzp,1:nxp);
[x2,z2] = meshgrid(1:nblock_p/nblock:nzp,1:nblock_p/nblock:nxp);
disp(['Actual size of histogram: ~',num2str(nblock_p),' pixels'])


%Segment photos into patches of histograms and then interpolate the actual
%histogram into the size of the simulated histograms(54 x 54 in this
%study).
dpixel = 2*npin*nblock_p/nblock; %how many pixels for one grid of translation. 

for pn = 1:repnum*circnum^2 %Process circnum x circum photos with "repnum" repetitions.
    
    bkg = double((imread([photopath,'/',list1(2*pn-1).name])));
    blur = double((imread([photopath,'/',list1(2*pn).name])));
    %load FUS photo(blur) and its background photo (bkg) for
    %segmentation.  
    %Note that here the transducer is on the right side of photo. 
    
    ss =order(mod(pn-1,circnum^2)+1)-1; % This is the ss-th photo in the current repetition. 
    xoffset = round(mod(ss,circnum)*dpixel); 
    zoffset = round(floor(ss/circnum)*dpixel);
    
    %Take out the green channel of RGB photo, move the rectangular FOV by xoffset pixels and zoffset pixels following the photos'
    %grid translations, and crop the photo into the smaller rectangular
    %photo.
    bkg = (bkg(xxmin+xoffset:xxmax+xoffset,zzmin+zoffset:zzmax+zoffset,2));
    blur =(blur(xxmin+xoffset:xxmax+xoffset,zzmin+zoffset:zzmax+zoffset,2));
    
    %Interpolate the actual histograms (nblock_p x nblock_p)to the patches of the same size as
    %simulated histograms (nblock x nblock)
    bkg = interp2(x1,z1,bkg,x2,z2,'linear');
    blur = interp2(x1,z1,blur,x2,z2,'linear');
    bkg(bkg<0) = 0;
    blur(blur<0) = 0;

    %Preprocess the non-FUS photo to have "bwmorph" function find the edges
    %of each histogram. 
    bkgtmp =  imguidedfilter(filter2(fspecial('average',3),(bkg/max(bkg(:)))));
    bkgtmp = (bkgtmp-min(bkgtmp(:)))/(max(bkgtmp(:))-min(bkgtmp(:)));
    [N,X] = hist(bkgtmp(:));
    edges = filter2(fspecial('average',5),(bwmorph(bkgtmp>median(X),'shrink',inf) == 1))>0.1;
    
    %Find the coordinates(xx,zz) of the center of histograms in the photo. 
    [m,n] = find(edges == 1); %#ok<IM2BW>
    tbl1 = tabulate(m);
    tbl2 = tabulate(n);
    thres1 = smoothdata(tbl1(:,2));
    thres2 = smoothdata(tbl2(:,2));
    thres1 = thres1/max(thres1(:));
    thres2 = thres2/max(thres2(:));
    [~,zz] = findpeaks(thres2,'MinPeakHeight',0.7);
    [~,xx] = findpeaks(thres1,'MinPeakHeight',0.7);
    nxx = length(xx);
    nzz = length(zz);
    zzind = find(diff(zz)<(nblock*0.95));
    xxind = find(diff(xx)<(nblock*0.95));
    zzind(zzind>mean(zzind)) = zzind(zzind>mean(zzind))+1;
    xxind(xxind>mean(xxind)) = xxind(xxind>mean(xxind))+1;
    xx(xxind) = [];
    zz(zzind) = [];
    xx = round((xx(1:end-1)+xx(2:end))/2);
    zz = round((zz(1:end-1)+zz(2:end))/2);
    xx(xx<(round(nblock/2)+1)) = [];
    xx(xx>(size(bkg,1)-round(nblock/2))) = [];
    zz(zz<(round(nblock/2)+1)) = [];
    zz(zz>(size(bkg,2)-round(nblock/2))) = [];
    
    %By the above process, we can usually get the same number of segmented
    %histograms for the actual photos with different grid translations. 
    display([num2str(ss),'-th,# of histograms in x-dimension:',num2str(length(xx)),',# of histograms in z-dimension:',num2str(length(zz))]);
    
    
    histblock = [];
    unblurblock = [];
    %Scale the green-channel photo (0~255) into 0~1.
    bkg = bkg/255;
    blur = blur/255;
    for ii = 1:length(zz)
        for jj = 1:length(xx)
            
            %Basically, the center coordinates calculated above are not
            %accurate. It is corrected by registration of non-FUS histograms with the very simple simulated
            %histogram given at the beginning of this script.
            tmpmin = imguidedfilter(imguidedfilter((bkgtmp(xx(jj)-nhalf+1:xx(jj)+nhalf,zz(ii)-nhalf+1:zz(ii)+1+nhalf))));
            tmpmin = (tmpmin-min(tmpmin(:)))/(max(tmpmin(:))-min(tmpmin(:)));
            tmpmin = tmpmin>graythresh(tmpmin);
            c = xcorr2(1-tmpmin,1-blocks);
            c = c(nhalf+1:end-nhalf,nhalf+1:end-nhalf);
            [mm,nn] = find(c == max(c(:)));
            minm = round(mean(mm))-nhalf;
            minn = round(mean(nn))-nhalf;
            tmpunblur = (bkg(xx(jj)+minm-nhalf+1:xx(jj)+minm+nhalf,zz(ii)+minn-nhalf+1:zz(ii)+minn+nhalf));
            tmpblur = (blur(xx(jj)+minm-nhalf+1:xx(jj)+minm+nhalf,zz(ii)+minn-nhalf+1:zz(ii)+minn+nhalf));
            unblurblock(:,:,jj,ii) = tmpunblur;
            histblock(:,:,jj,ii) =  tmpblur;
        end
    end
    nX = size(histblock,3);
    nZ = size(histblock,4);
    seg.blur = histblock;
    seg.unblur = unblurblock;
    seg.nx = nX;
    seg.nz = nZ;
    save([datapath,'/seg_',num2str(ss+1),'_',num2str(floor((pn-1)/length(order))+1),'.mat'],'seg','xxmin','xxmax','zzmin','zzmax');
end


%%Then the segmented histograms are tiled into the interleaved photo
%%achieving high resolution. "tilehisblur" is used to reconstruct RMS
%%projected pressure. 
tilehisblur = [];
tilehisunblur = [];
for  ff = 1:repnum
    for nn = 1:circnum^2
        load([datapath,'/seg_',num2str(nn),'_',num2str(ff),'.mat']);
        nff = mod(nn-1,circnum)+1;
        mff = (floor(((nn)-1)/circnum)+1);
        nX = seg.nx;
        nZ = seg.nz;
        disp([num2str(nff),',',num2str(mff),',',num2str(nX),',',num2str(nZ)])
        tilehisblur(:,:,nff:circnum:nX*circnum,mff:circnum:nZ*circnum) = seg.blur;
        tilehisunblur(:,:,nff:circnum:nX*circnum,mff:circnum:nZ*circnum) = seg.unblur;
    end
    nX = size(tilehisblur,3);
    nZ = size(tilehisblur,4);
    tilehisblur = single(tilehisblur);
    tilehisunblur = single(tilehisunblur);
    save([datapath,'/tilehis',num2str(ff)],'tilehisblur','tilehisunblur','nX','nZ'); % tilehis*.mat will finally be brought into the input of the neural network.
end



