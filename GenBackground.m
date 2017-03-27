%% 
clear
clc
close all;

piIP = '10.114.29.211'; 
piPW = 'zebrazebra'; % password for the pi
prefix = 'image30'; 
capturecmd = ['fswebcam -d /dev/video0 -r 2304x1536 -p YUYV -s brightness=128 --no-banner --png 2 ' ...
    prefix '.png'];
combinedcmd = ['/usr/local/bin/sshpass -p ' piPW ' ssh pi@' piIP ' ' capturecmd];
dlcmd = ['/usr/local/bin/sshpass -p ' piPW ' scp pi@' piIP ':~/' prefix '.png ./'];

bgimg='img';
ulcmd = ['/usr/local/bin/sshpass -p ' piPW ' scp /Users/masonchen/Downloads/Zebra_Code/' bgimg '.png pi@' piIP ':~/'];
displaycmd = ['sudo fbi -d /dev/fb0 -T 1 -noverbose ' bgimg '.png'];
combcmd = ['/usr/local/bin/sshpass -p ' piPW ' ssh pi@' piIP ' ' displaycmd];


%% Zebra

    wavNoiseImg = ones(1200,1920);
    lineGap = 20; % pixels between lines
            nZeros = 2;
            template = ones(lineGap,1);
            template(1:nZeros) = 0;
            cols = kron(ones(floor(size(wavNoiseImg,2)/lineGap),1),template);
            cols = [cols;ones(size(wavNoiseImg,2)-length(cols),1)]; 
            wavNoiseImg = ones(size(wavNoiseImg,1),1)*cols.';
            img=wavNoiseImg;
            system(ulmd);
            system(combcmd);
            
            system(combinedcmd);
            system(dlcmd);
            offpic = imread([prefix '.png']);
            picoff=0;
%             for reps=1:total_reps %for-loop for total_reps to take that many # of pictures
%                 imge1 = snapshot(cam);
%                 off_matrix(:,:,reps)=imge1(:,:,1);
%                 picoff=double(rgb2gray(imge1))+double(picoff);
%             end
%             picoff=picoff./total_reps;
%             picoff=picoff-min(picoff(:));
%             picoff=picoff/max(picoff(:));
            off = 'off';
%             if length(angles_vector)>1
%                str2 = [str '/' angle '.' off '.jpg'];
%             else
%                str2 = [str '/' off '.jpg'];
%             end
%             imwrite(picoff,str2);
            picon=0;
%             for reps=1:total_reps
%                 imge1 = snapshot(cam);
%                 on_matrix(:,:,reps)=imge1(:,:,1);
%                 picon=double(rgb2gray(imge1))+double(picon);
%             end
%             picon=picon./total_reps;
%             picon=picon-min(picon(:));
%             picon=picon/max(picon(:));
            on = 'on';
            system(combinedcmd);
            system(dlcmd);
            onpic = imread([prefix '.png']);
%             if length(angles_vector)>1
%                 str3 = [str '/' angle '.' on '.jpg'];
%             else
%                 str3 = [str '/' on '.jpg'];
%             end;
%             imwrite(picon,str3);
            sub=onpic-offpic;
%             if length(angles_vector)>1
%                 str4 = [str '/' angle '.sub' '.jpg'];
%             else
%                 str4 = [str '/' 'sub' '.jpg'];
%             end;
%             sub=sub-min(sub(:));
%             sub=sub/max(sub(:));
%             imwrite(sub,str4);
 

%% Hex color

    img = ones(1200,1920);
    %%% change rectangle size below %%%
     width=3; %pixels
     leng=3; %pixels
     gap=12;
     rows=floor(1200/gap);
     cols=floor(1920/gap);
    for ii=0:2:rows-1
        for jj=0:2:cols-1
            img(ii*gap+1:ii*gap+leng,jj*gap+1:jj*gap+width)=2;
        end
        for jj=1:2:cols-1
             img(ii*gap+1:ii*gap+leng,jj*gap+1:jj*gap+width)=4;
        end
     end
     for ii=1:2:rows-1
         for jj=1:2:cols-1
            img(ii*gap+1:ii*gap+leng,jj*gap+1:jj*gap+width)=2;
         end
         for jj=0:2:cols-1
        img(ii*gap+1:ii*gap+leng,jj*gap+1:jj*gap+width)=3; 
         end
     end

       map=[0,0,0;0,0,1;0,1,0;1,0,0];
        
%             system(ulcmd);
%             system(combcmd);

%% Hex
     img = ones(1200,1920);
    %%% change rectangle size below %%%
     width=3; %pixels
     leng=3; %pixels
     gap=12;
     rows=floor(1200/gap);
     cols=floor(1920/gap);
    for i=0:2:rows-1
        for j=0:2:cols-1
            img(i*gap+1:i*gap+leng,j*gap+1:j*gap+width)=0;
        end
    end
        for ii=1:2:rows-1
            for jj=1:2:cols-1
        img(ii*gap+1:ii*gap+leng,jj*gap+1:jj*gap+width)=0;
            end
        end
  figure;
  imagesc(img); colormap gray
        hold on;
        axis off; %turn off axes of figure
        
        k=1;
        bgimg{k}=img;
        for colshift=0:width:gap-width
               for rowshift=0:leng:2*gap-leng
                bgimg{k}=imtranslate(img,[colshift,rowshift]);
                imagesc(bgimg{k}); colormap(gray);
                R=sprintf('hex%d',k);
                imwrite(bgimg{k},[R, '.png']);
                k=k+1;
               end
        end
       
