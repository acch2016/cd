clear all; clc;
% diseño del filtro
m = [1 1 0 0];
o = 100;
f = [0 0.05 0.05 1]
f1 = fir2(o,f,m);
% fvtool(f1)
f = [0 0.2 0.2 1]
f2 = fir2(o,f,m);
% fvtool(f2)
f = [0 0.4 0.4 1]
f3 = fir2(o,f,m);
% fvtool(f3)
%% filter
% ventajas vectores mismo tamaño
% xnrz_filtered = filter(f1,1,xnrz);% denominador 1
%% conv
% en la convolucion suma de los dos vectores -1
% si deseo observar transitorios utilizo convolucion
% xnrz_filtered = conv(f1,xnrz)
%%  Unipolar NRZ % % %
bits = randi([0 1],1024,1); % Generarte random bits
mp = 10; % samples per pulse
Fs = 96000;
Ts = 1/Fs;
% The bit rate is Rb= Rs= Fs / mp, because 1 bit= 1 %symbol and % every symbol has mpsamples per bit
Rs = Fs / mp;
pbase = rectwin(mp); % Pulso completo;
% wvtool(pbase)
s = zeros(1,numel(bits)*mp);
s(1:mp:end) = bits;
% subplot(3,1,1); stem(s(1:mp*16)); title('tren de impulsos');
unrz = conv(pbase,s);
% subplot(3,1,3);pwelch(unrz,[],[],[],Fs,'power'); title('PSD of Unipolar NRZ');
% subplot(3,1,2);plot(unrz(1:mp*16)); title('Unipolar NRZ');
%% Polar NRZ % 
mp = 10; % samples per pulse
Fs = 96000;
Ts = 1/Fs;
% The bit rate is Rb= Rs= Fs / mp, because 1 bit= 1 %symbol and % every symbol has mpsamples per bit
Rs = Fs / mp;
pbase = rectwin(mp); % Pulso completo;
% wvtool(pbase);
sym = (bits * 2)-1;
s = zeros(1,numel(bits)*mp);
s(1:mp:end) = sym;
% subplot(3,1,1); stem(s(1:mp*16)); title('tren de impulsos');
pnrz = conv(pbase,s);
% subplot(3,1,3);pwelch(pnrz,[],[],[],Fs,'power'); title('PSD of Polar NRZ');
% subplot(3,1,2);plot(pnrz(1:mp*16)); title('Polar NRZ');
%% Bipolar NRZ % 
mp = 10; % samples per pulse
Fs = 96000;
Ts = 1/Fs;
% The bit rate is Rb= Rs= Fs / mp, because 1 bit= 1 %symbol and % every symbol has mpsamples per bit
Rs = Fs / mp;
pbase = rectwin(mp); % Pulso completo;
wvtool(pbase)
am = mod(1:length(bits(bits == 1)), 2); % 0s y 1s del tamaño de los bits en 1 en la variable bits
am(am == 0) = -1; % alternate mark (vector de -1s y 1s)
bits1 = bits; % para no modificar la variable original que contine los bits
bits1(bits == 1) = am; % alternate mark inversion AMI (tres estados)
s = zeros(1,numel(bits1)*mp);
s(1:mp:end) = bits1;
subplot(3,1,1);stem(s(1:mp*16)); title('tren de impulsos');
bnrz = conv(pbase,s);
subplot(3,1,3);pwelch(bnrz,[],[],[],Fs,'power'); title('PSD of Bipolar NRZ');
subplot(3,1,2);plot(bnrz(1:mp*16));title('Bipolar NRZ');
figure;powerbw(bnrz,Fs) % para obtener el bw de -3dB
%% Manchester % %
bits = randi([0 1],1024,1); % Generarte random bits
mp = 10; % samples per pulse
Fs = 96000;
Ts = 1/Fs;
% The bit rate is Rb= Rs= Fs / mp, because 1 bit= 1 %symbol and % every symbol has mpsamples per bit
Rs = Fs / mp
Rb = Rs
pM = [-ones(1,mp/2) ones(1,mp/2)]; %pulso base
% wvtool(pM);
sym = (bits * 2)-1;
s = zeros(1,numel(bits)*mp);% tren de impulsos
s(1:mp:end) = sym; %tren
% subplot(3,1,1); stem(s(1:mp*16)); title('tren de impulsos'); 
xM = conv(pM,s);
% subplot(3,1,3); pwelch(xM,200,100,200,Fs,'power'); title('PSD of Manchester');
% subplot(3,1,3); pwelch(xM,[],[],[],Fs,'power'); title('PSD of Manchester');
% subplot(3,1,2); plot(xM(1:mp*16));title('Manchester');
%%
unrz_filtered = filter(f1,1,unrz);% denominador 1
subplot()
pwelch(unrz,[],[],[],Fs,'power'); title('PSD of filtered Unipolar NRZ');