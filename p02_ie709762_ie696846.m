%% PARTE I
clear;clc
%% Tx
Fs = 192e3;
Ts = 1/Fs;
t = 0:Ts:1-Ts;
signal = sin(2*pi*5000*t); 
soundsc(signal, Fs);
%% Rx
info = audiodevinfo
info.input.ID
recorder1 = audiorecorder(192000,16,1,1)
record(recorder1);
stop(recorder1)
play(recorder1)
audioarray = getaudiodata(recorder1);
t = 0:1/192000:(numel(audioarray)-1)/192000;
plot(audioarray)
%% Impulso
y = [zeros(1,Fs) 1 zeros(1,Fs)];
t = 0:Ts:2;
plot(t,y)
soundsc(y,192000);
wvtool(y)
%% AWGN
B = 1e6;
T = 1;
Ts = 1 / ( 2 *B ) ;
t = 0 : Ts : T-Ts ;
M = 2*B*T;
% y = sqrt ( potr )* randn ( 1 ,M ) ;
y = randn ( 1, M ) ;
% pn = sum(y.^2)/numel(y) 
% y = y/sqrt(pn);
% pn = sum(y.^2)/numel(y)
soundsc(y,Fs);
% figure('Position',[100, 100, 250, 250]); plot(t,y)
figure('Position',[100, 100, 500, 200]); pwelch(y,500,300,500,'one-side', 'power', Fs)
%% Chirp
clc;clear
Fs = 44.1e3;
Ts =1/Fs;
t = 0:Ts:5;
y = chirp(t,100,5,20000,'linear');
soundsc(y,Fs);
figure('Position',[100, 100, 500, 200]); pspectrum(y,Fs,'spectrogram','TimeResolution',0.2, 'OverlapPercent',99,'Leakage',0.85)
figure('Position',[100, 100, 500, 200]); pwelch(y,500,300,500,'one-side', 'power', Fs)
%% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % FASE II Tx
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
clc; clear;
%% Diseño del pulso base SRRC
beta = 0.5;
B = 1.8e3;
Rb = 2*B /(1+beta);
fs = 192000;
% fs = 96000;
mp = fs/Rb;
% Rb = fs/mp; % bps
D = 6;
Tp=1/Rb;
Ts=1/fs;
type = 'srrc';
energy = Tp;
[p t] = rcpulse(beta,D,Tp,Ts,type,energy);
% wvtool(p)
% fvtool(p)
% plot(t,p)
e = Ts*p*p'
%% LENA 
load lena512.mat
% whos('-file','lena512.mat')
% imshow(uint8(lena512))
lenarec=lena512(252:284,318:350);
% imshow(uint8(lenarec))
b=de2bi(lenarec,8);
% b=de2bi(lenarec,8,'left-msb'); % % % % % % % % % % % % % version profesor
b=b';
bits=b(:);
%% [ Preambulo Header Payload ]
pr = [1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 1 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 1];
altoancho =  de2bi(size(lenarec), 8, 'left-msb');
header = altoancho(1,:);
header(1,9:16) = altoancho(2,:);
bits = [pr'; header'; bits];
%% tren de pulsos pnrz
sym = bits*2-1; % symbols
s = zeros(1,numel(bits)*mp);
s(1:mp:end) = sym; % tren de impulsos
pnrz = conv(p, s); %tren de pulsos
% figure('Position',[100, 100, 700, 250]); plot(pnrz(1:10*mp));
% title('Transmitted signal')
% ylabel('Amplitude (V)')
% xlabel('Time (s)')
% plot(pnrz);
%% Diagrama de ojo
eyediagram(pnrz(3*mp:603*mp),3*mp)
%% PSD
figure('Position',[100, 100, 700, 250]); pwelch(pnrz(32+16:end),500,300,500,'one-side', 'power', fs);
%% Tx
soundsc([zeros(1, fs/2) pnrz zeros(1, fs/2)],fs);
%% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % FASE III Tx
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
clc; clear
%% Diseño del pulso base SRRC
beta = 0.5;
B = 1e3;
Rb = 2*B /(1+beta);
fs = 96000;
% fs = 96000;
mp = fs/Rb;
% Rb = fs/mp; % bps
D = 6;
Tp=1/Rb;
Ts=1/fs;
type = 'srrc';
energy = Tp;
[p t] = rcpulse(beta,D,Tp,Ts,type,energy);
wvtool(p)
% fvtool(p)
% plot(t,p)
e = Ts*p*p'
%% LENA 
load lena512.mat
% whos('-file','lena512.mat')
% imshow(uint8(lena512))
lenarec=lena512(252:284,318:350);
% imshow(uint8(lenarec))
b=de2bi(lenarec,8);
% b=de2bi(lenarec,8,'left-msb'); % % % % % % % % % % % % % version profesor
b=b';
bits=b(:);
%% [ Preambulo Header Payload ]
pr = [1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 1 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 1];
altoancho =  de2bi(size(lenarec), 8, 'left-msb');
header = altoancho(1,:);
header(1,9:16) = altoancho(2,:);
bits = [pr'; header'; bits];
%% tren de pulsos pnrz
sym = bits*2-1; % symbols
s = zeros(1,numel(bits)*mp);
s(1:mp:end) = sym; % tren de impulsos
pnrz = conv(p, s); %tren de pulsos
% figure('Position',[100, 100, 700, 250]); plot(pnrz(1:10*mp));
% title('Transmitted signal')
% ylabel('Amplitude (V)')
% xlabel('Time (s)')
% plot(pnrz);
%% AM
Fc = 6.5e3; Fs = 96e3;
pnrz_am = ammod(pnrz,Fc,Fs); % Modulate.
% figure('Position',[100, 100, 700, 250]); pwelch(pnrz,500,300,500,'one-side', 'power', fs);
% figure('Position',[100, 100, 700, 250]); pwelch(pnrz_am,500,300,500,'one-side', 'power', Fs);
% fvtool([num,den]);
%% Tx
soundsc([zeros(1, fs/2) pnrz_am zeros(1, fs/2)], fs);