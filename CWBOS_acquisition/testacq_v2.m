%|Script to control BOS system to acquire FUS-photos with remotely
%|controlling waveform generator. 
%|
%|User can adjust parameters of waveform generator in this script. 



GenAddress='10.114.17.161'; % IP address of waveform generator
hornSchunck=false; %whether or not to execute HS method
%prevent from sending too much power to Xducer

%% TCP/IP connection for function generator
% Create TCP/IP object 'fncngen'. Specify server machine and port number.
fncngen = tcpip(GenAddress, 5025,'NetworkRole','Client');
% Set size of receiving buffer, if needed.
set(fncngen, 'InputBufferSize', 30000);

disp('opening connection..')
% Open connection to the server.
fopen(fncngen);
disp('connection created!');

% Initialize pulser
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
freq = 1.16; % MHz
amplitude = 200; %mVpp
if amplitude > 400
    error('amplitude is too high - 500mV is maximum')
end
fprintf(fncngen,'OUTP1 OFF;');
fprintf(fncngen,'OUTP1:LOAD 50.0');
fprintf(fncngen,'OUTP1:POL NORM');
mVpp=amplitude/1000;
mVpp2 = mVpp/2;
mVpp3 = -mVpp2;
fprintf(fncngen,'SOUR1:FUNC:SHAP SIN;');
cmd = sprintf('SOUR1:FREQ %1.3fe+06;',freq);
cmd2 = sprintf('SOUR1:VOLT %4.3f;',mVpp);
cmd3 = sprintf('SOUR1:VOLT:HIGH %4.3f;',mVpp2);
cmd4 = sprintf('SOUR1:VOLT:LOW %4.3f;',mVpp3);
fprintf(fncngen,cmd);
fprintf(fncngen,cmd2);
fprintf(fncngen,cmd3);
fprintf(fncngen,cmd4);
fprintf(fncngen,'SOUR1:VOLT:UNIT VPP;');
fprintf(fncngen,'SOUR1:VOLT:OFFS 0.0E+00;');
fprintf(fncngen,'UNIT:ANGLe DEG;');
fprintf(fncngen,'SOUR1:SWEep:STATe OFF;');
fprintf(fncngen,'SOUR1:SWEep:SPAC LIN;');
fprintf(fncngen,'SOUR1:SWEep:RTIMe 0.0E+00;');
fprintf(fncngen,'SOUR1:SWEep:HTIMe 0.0E+00;');
fprintf(fncngen,'SOUR1:FREQ:STOP 1.0E+03;');
fprintf(fncngen,'SOUR1:FREQ:STAR 1.0E+02;');
fprintf(fncngen,'UNIT:ANGLe DEG;');
fprintf(fncngen,'OUTP1 OFF;');
fprintf(fncngen,'OUTP1:LOAD 50.0');
fprintf(fncngen,'OUTP1:POL NORM');


%% Take pictures!
num_repetitions = 5; %number of repetitions of "on" (and "off") images that are averaged
a = serial('COM8','BaudRate',9600);
ip_address = '10.248.218.47'; %iPad address shown on the ipad screen in the app. 
dd = [0 2 4 6 16 18 20 22 32 34 36 38 48 50 52 54];
dd = reshape(dd, 4,4)';
for jj = 1:num_repetitions
    %Repeat acquistions 5 times.
    for ii =  1:length(dd(:))
        system(['py BOSTomoController.py p 1 8 ',num2str(dd(ii)),' ',num2str(0),' ',ip_address]);
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
%%
fprintf(fncngen,'OUTP1 OFF;');
