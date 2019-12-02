%|Script to control BOS system to acquire FUS-photos without remotely
%|controlling waveform generator. 
%|User needs to adjust parameters of waveform generator mannualy. 


delete(instrfindall); 
% Close all connections and deletes the objects at first.
% Need to adjust settings of waveform generator mannually. 
  
a = serial('COM8','BaudRate',9600);
ip_address = '10.248.218.47'; %iPad address shown on the ipad screen in the app. 
dd = [0 2 4 6 16 18 20 22 32 34 36 38 48 50 52 54];
dd = reshape(dd, 4,4)';
for jj = 1:5 
    %Repeat acquistions 5 times.
    for ii =  1:length(dd(:))
        system(['py BOSTomoController.py p 1 8 ',num2str(dd(ii)),' ',num2str(0),' ', ip_address]);
        pause(2);
        fopen(a);
        fprintf(a,'%i',0);
        fprintf(a,'%i',0);
        fclose(a);
        pause(2);  %Wait for the camera uploads the taken photos to the experment computer. 
        fopen(a);
        fprintf(a,'%i',1);
        fprintf(a,'%f',50); 
        % With FUS on, send a TTL pulse to waveform generator, and wait 50ms and open the camera shutter. 
        fclose(a);
        pause(2); %Wait for the camera uploads the taken photos to the experment computer. 
    end
end

