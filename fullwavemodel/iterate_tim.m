function [recon,projnew,projprenew,diff,corr,acor] = iterate_tim(projori,projpre)

for i = 1:size(projori,3)
    for j = 1:size(projori,1)
        A = squeeze(projpre(j,:,i));
        B = squeeze(projori(j,:,i));
        xx = corrcoef(A,B);    
        corr(i,j) = xx(2);
        if(isnan(xx(2)))
        corr(i,j) = 0;
        end
        
  
        [acor(i,j,:),lag] = xcorr(B,A);
        [~,I] = max(abs(squeeze(acor(i,j,:))));
%         d = round(lag(I)/2);
        d = round(KL(A,B,5,5,1));
        if (isnan(d)||d<0)
            d = 0;
        end
        tmp = [B,B];
        tmp = real(ifft(fft(tmp).*exp(1i*d*(0:(length(tmp)-1))*2*pi/length(tmp))));
        diff(i,j) = d;
        projnew(j,:,i) = tmp(1:length(B)); 

%         if i == 71 && j==2
%         display(d);
%         plot(squeeze(projnew(j,:,i)),'--')
%         hold on;
%         pause(0.000005)
%         end
    end
end
for i = 1:size(projori,2)
    recon(:,:,i)= recon_eachpoint(squeeze(projnew(:,i,:)));
    % figure; im(drk);
    % figure; im(P(:,:,153,12));
end

P = recon;
angles = 0:179;
kk = -1:2/(size(P,1)-1):1;
kk = [kron(cosd(angles(:)),kk.') kron(sind(angles(:)),kk.')];
G = Gmri(kk,true(size(P,1)),'fov',size(P,1)/2);
for j = 1:size(P,3)
    I = reshape(G * (squeeze(P(:,:,j))),[size(P,1),size(angles)]);
    projprenew(:,j,:) = squeeze(real(ifftshift(ifft(ifftshift(I,1),[],1))))';   
end


