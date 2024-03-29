%% 
bits = randi([0 1],1024,1); % Generarte random bits
mp = 20; % samples per pulse
Fs = mp; % Can be different
Ts=1/Fs;
% The bit rate is Rb= Rs= Fs / mp, because 1 bit= 1 %symbol and % every symbol has mpsamples per bit
% unipolar nrz
%pnrz = [1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1]; % Pulse type with mp samples
pnrz = 5* ones(1,mp); % Similar to Arduino pulses with 5V % pulso base
s = zeros(1,numel(bits)*mp);
s(1:mp:end) = bits;
stem(s(1:mp*16))
wvtool(pnrz)
xnrz = conv(pnrz,s);
figure;pwelch(xnrz,[],[],[],Fs,'power'); %title('PSD of Unipolar NRZ');
figure;plot(xnrz(1:mp*16))
%% prz
clc;clear all;
bits = randi([0 1],1024,1); % Generarte random bits
mp = 20; % samples per pulse
Fs = mp; % Can be different
Ts=1/Fs;
% The bit rate is Rb= Rs= Fs / mp, because 1 bit= 1 %symbol and % every symbol has mp samples per bit
% polar rz
pRZ = 5 * [ones(1,mp/2) zeros(1,mp/2)]
sym = (bits * 2)-1
s = zeros(1,numel(bits)*mp);
% s(1:mp:end) = bits;
s(1:mp:end) = sym;
stem(s(1:mp*16))
wvtool(pRZ)
xRZ = conv(pRZ,s);
figure;pwelch(xRZ,[],[],[],Fs,'power'); %title('PSD of Unipolar NRZ');
figure;plot(xRZ(1:mp*16))
%% Manchester
bits = randi([0 1],1024,1); % Generarte random bits
mp = 20; % samples per pulse
Fs = mp; % Can be different
Ts=1/Fs;
% The bit rate is Rb= Rs= Fs / mp, because 1 bit= 1 %symbol and % every symbol has mp samples per bit
% Manchester
pM = 5 * [-ones(1,mp/2) ones(1,mp/2)] %pulso base
sym = (bits * 2)-1
s = zeros(1,numel(bits)*mp);% tren de impulsos
% s(1:mp:end) = bits;
s(1:mp:end) = sym;%tren
stem(s(1:mp*16))
wvtool(pM)
xM = conv(pM,s);
figure;pwelch(xM,[],[],[],Fs,'power'); %title('PSD of Unipolar NRZ');
 plot(xM(1:mp*32))