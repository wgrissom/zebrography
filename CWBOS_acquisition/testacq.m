%|Matlab Script to control BOS system to acquire FUS-photos without remotely
%|controlling waveform generator. 
%|User needs to adjust parameters of waveform generator mannualy. 


delete(instrfindall); 
% Close all connections and deletes the objects at first.
% Need to adjust settings of waveform generator mannually. 
  
a = serial('COM8','BaudRate',9600); % creat the serial commucation between Matlab and Arduino. 
ip_address = '10.248.218.47'; % iPad address shown on the ipad screen in the app. 
dd = [0 2 4 6 16 18 20 22 32 34 36 38 48 50 52 54]; % offset of background pattern. Here, pins are shifted by 2n pixels at each step down first and then right. 
dd = reshape(dd, 4,4)';

for jj = 1:5 
    %Repeat acquisitions 5 times for each background pattern
    for ii =  1:length(dd(:))
        system(['py BOSTomoController.py p 1 8 ',num2str(dd(ii)),' ',num2str(0),' ', ip_address]); % send the command to iPad
        % p represents "pins"; pin width is 1 pixel and pin spacing is 8 pixels; 
        % 0 represents white background (1: black background); ip_address: IP address displayed on the iPad
        pause(2); 
        fopen(a); % open the serial communication
        fprintf(a,'%i',0); % output of waveform generator is off. 
        fprintf(a,'%i',0); % tell the Arduino Board when to trigger the waveform generator. 
        fclose(a);
        pause(2);  % wait for the camera to upload the taken photos to the experiment computer. Not sure if the serial communication needs to be turned off and on again here. 
        fopen(a);
        fprintf(a,'%i',1); % trigger the waveform generator to send a pulse
        fprintf(a,'%f',50); % tell the Arduino Board when (wait for 50 ms) to trigger the waveform generator. 
        % With FUS on, open the camera shutter, wait 50ms and send a TTL pulse to waveform generator. 
        fclose(a); % close the serial communication
        pause(2); % wait for the camera to upload the taken photos to the experiment computer. 
    end
end

