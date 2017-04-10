close all; clc;
clear all
IPaddress='10.114.29.34'; 
amplitude=200; 
total_reps=10; %number of repetitions of "on" (and "off") images that are averaged
hornSchunck=false; %whether or not to execute HS method

%prevent from sending too much power to Xducer
if amplitude > 400
    error('amplitude is too high - 500mV is maximum')
end
        
%make strings of time for purposing of naming saved files chronologically
time=clock;
year = num2str(time(1));
month = num2str(time(2));
day = num2str(time(3));

%If there is no data storage location in current folder, create one
str = [pwd '\Data'];
if exist(str) ~= 7
    mkdir(str)
end 
       
%If there is not a data storage location for today's data, create one
str = [str '/' year '.' month '.' day];
if exist(str) ~= 7
    mkdir(str)
end

%% TCP/IP connection
% Create TCP/IP object 'fncngen'. Specify server machine and port number. 
fncngen = tcpip(IPaddress, 5025,'NetworkRole','Client'); 

% Set size of receiving buffer, if needed. 
set(fncngen, 'InputBufferSize', 30000); 

disp('opening connection..')
% Open connection to the server. 
fopen(fncngen);
disp('connection created!');

% Initialize pulser
freq = 1.16; % MHz 
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

%% Create server (run by "evaluate selection")
    import java.net.ServerSocket
    import java.io.*

    server_socket  = [];
    output_socket  = [];
    output_port= input('Enter port numer: '); %%set port here
    message=' connected'; %write message here
    server_socket = ServerSocket(output_port);
    output_socket = server_socket.accept;
    fprintf(1, 'Client connected\n');
    
    output_stream   = output_socket.getOutputStream;
    d_output_stream = DataOutputStream(output_stream);
    fprintf(1, 'Writing %d bytes\n', length(message));
    d_output_stream.writeBytes(char(message));
    %pause(3);
    
%% 
new_message=' gap12';
fprintf(1, 'Writing %d bytes\n', length(new_message));
d_output_stream.writeBytes(char(new_message));
d_output_stream.flush;
fprintf(fncngen,'OUTP1 OFF;');
%TAKE A PICTURE!!!

%% 
import java.awt.Robot;
import java.awt.event.*;
mouse = Robot;
screenSize = get(0, 'screensize');

for ii=1:32
    new_message=sprintf(' hex%d',ii);
    fprintf(1, 'Writing %d bytes\n', length(new_message));
    d_output_stream.writeBytes(char(new_message));
    pause(5);

    mouse.mousePress(InputEvent.BUTTON1_MASK);
    mouse.mouseRelease(InputEvent.BUTTON1_MASK);
    pause(2);
    
    fprintf(fncngen,'OUTP1 ON;');
    mouse.mousePress(InputEvent.BUTTON1_MASK);
    mouse.mouseRelease(InputEvent.BUTTON1_MASK);
    pause(0.5);

    fprintf(fncngen,'OUTP1 OFF;');
end
%% Close port
server_socket.close;
output_socket.close;

%% Rename pictures
% Get all PDF files in the current folder
files = dir('*.pdf');
% Loop through each
for id = 1:length(files)
    % Get the file name (minus the extension)
    [~, f] = fileparts(files(id).name);
      % Convert to number
      num = str2double(f);
      if ~isnan(num)
          % If numeric, rename
          movefile(files(id).name, sprintf('%03d.pdf', num));
      end
end


