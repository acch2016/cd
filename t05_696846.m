clear all; clc;
% diseño del filtro
m = [1 1 0 0];
o = 100;
f = [0 0.05 0.05 1];
f1 = fir2(o,f,m);
% fvtool(f1)
f = [0 0.2 0.2 1];
f2 = fir2(o,f,m);
% fvtool(f2)
f = [0 0.4 0.4 1];
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
% subplot(3,1,3);
figure;pwelch(unrz,[],[],[],Fs,'power'); title('PSD of Unipolar NRZ');
% subplot(3,1,2);plot(unrz(1:mp*16)); title('Unipolar NRZ');
%% Polar NRZ % 
% wvtool(pbase);
sym = (bits * 2)-1;
s = zeros(1,numel(bits)*mp);
s(1:mp:end) = sym;
% subplot(3,1,1); stem(s(1:mp*16)); title('tren de impulsos');
pnrz = conv(pbase,s);
% subplot(3,1,3);
figure;pwelch(pnrz,[],[],[],Fs,'power'); title('PSD of Polar NRZ');
% subplot(3,1,2);plot(pnrz(1:mp*16)); title('Polar NRZ');
%% Bipolar NRZ % 
% wvtool(pbase)
am = mod(1:length(bits(bits == 1)), 2); % 0s y 1s del tamaño de los bits en 1 en la variable bits
am(am == 0) = -1; % alternate mark (vector de -1s y 1s)
bits1 = bits; % para no modificar la variable original que contine los bits
bits1(bits == 1) = am; % alternate mark inversion AMI (tres estados)
s = zeros(1,numel(bits1)*mp);
s(1:mp:end) = bits1;
% subplot(3,1,1);stem(s(1:mp*16)); title('tren de impulsos');
bnrz = conv(pbase,s);
% subplot(3,1,3);
figure;pwelch(bnrz,[],[],[],Fs,'power'); title('PSD of Bipolar NRZ');
% subplot(3,1,2);plot(bnrz(1:mp*16));title('Bipolar NRZ');
% figure;powerbw(bnrz,Fs) % para obtener el bw de -3dB
%% Manchester % %
pM = [-ones(1,mp/2) ones(1,mp/2)]; %pulso base
% wvtool(pM);
sym = (bits * 2)-1;
s = zeros(1,numel(bits)*mp);% tren de impulsos
s(1:mp:end) = sym; %tren
% subplot(3,1,1); stem(s(1:mp*16)); title('tren de impulsos'); 
xM = conv(pM,s);
% subplot(3,1,3); pwelch(xM,200,100,200,Fs,'power'); title('PSD of Manchester');
% subplot(3,1,3); 
figure;pwelch(xM,[],[],[],Fs,'power'); title('PSD of Manchester');
% subplot(3,1,2); plot(xM(1:mp*16));title('Manchester');
%% filtrar y graficar espectro
unrz_filter005 = filter(f1,1,unrz);% denominador 1
subplot(6,2,1)
pwelch(unrz_filter005,[],[],[],Fs,'power'); title('PSD Unipolar NRZ fc = 0.05');

pnrz_filter005 = filter(f1,1,pnrz);% denominador 1
subplot(6,2,2)
pwelch(pnrz_filter005,[],[],[],Fs,'power'); title('PSD Polar NRZ fc = 0.05');

bnrz_filter005 = filter(f1,1,bnrz);% denominador 1
subplot(6,2,3)
pwelch(bnrz_filter005,[],[],[],Fs,'power'); title('PSD Bipolar NRZ fc = 0.05');

xM_filter005 = filter(f1,1,xM);% denominador 1
subplot(6,2,4)
pwelch(xM_filter005,[],[],[],Fs,'power'); title('PSD Manchester fc = 0.05');


unrz_filter02 = filter(f2,1,unrz);% denominador 1
subplot(6,2,5)
pwelch(unrz_filter02,[],[],[],Fs,'power'); title('PSD Unipolar NRZ fc = 0.2');

pnrz_filter02 = filter(f2,1,pnrz);% denominador 1
subplot(6,2,6)
pwelch(pnrz_filter02,[],[],[],Fs,'power'); title('PSD Polar NRZ fc = 0.2');

bnrz_filter02 = filter(f2,1,bnrz);% denominador 1
subplot(6,2,7)
pwelch(bnrz_filter02,[],[],[],Fs,'power'); title('PSD Bipolar NRZ fc = 0.2');

xM_filter02 = filter(f2,1,xM);% denominador 1
subplot(6,2,8)
pwelch(xM_filter02,[],[],[],Fs,'power'); title('PSD Manchester fc = 0.2');


unrz_filter04 = filter(f3,1,unrz);% denominador 1
subplot(6,2,9)
pwelch(unrz_filter04,[],[],[],Fs,'power'); title('PSD Unipolar NRZ fc = 0.4');

pnrz_filter04 = filter(f3,1,pnrz);% denominador 1
subplot(6,2,10)
pwelch(pnrz_filter04,[],[],[],Fs,'power'); title('PSD Polar NRZ fc = 0.4');

bnrz_filter04 = filter(f3,1,bnrz);% denominador 1
subplot(6,2,11)
pwelch(bnrz_filter04,[],[],[],Fs,'power'); title('PSD Bipolar NRZ fc = 0.4');

xM_filter04 = filter(f3,1,xM);% denominador 1
subplot(6,2,12)
pwelch(xM_filter04,[],[],[],Fs,'power'); title('PSD Manchester fc = 0.4');
%% Grafique (utilizando la función stem ( ) ) las señales originales y las filtradas (16 figuras) usando solamente los 
% primeros 16 bits ( es decir con 16 x mp muestras + el retardo del filtro). 

figure;
subplot(4,1,1); stem(unrz(1:mp*16)); title('Unipolar NRZ');
subplot(4,1,2); stem(unrz_filter005(47:mp*16+47)); title('Unipolar NRZ fc = 0.05');
subplot(4,1,3); stem(unrz_filter02(50:mp*16+50)); title('Unipolar NRZ fc = 0.2');
subplot(4,1,4); stem(unrz_filter04(50:mp*16+50)); title('Unipolar NRZ fc = 0.4');
figure;
subplot(4,1,1); stem(pnrz(1:mp*16)); title('Polar NRZ');
subplot(4,1,2); stem(pnrz_filter005(47:mp*16+47)); title('Polar NRZ fc = 0.05');
subplot(4,1,3); stem(pnrz_filter02(50:mp*16+50)); title('Polar NRZ fc = 0.2');
subplot(4,1,4); stem(pnrz_filter04(50:mp*16+50)); title('Polar NRZ fc = 0.4');
figure;
subplot(4,1,1); stem(bnrz(1:mp*16));title('Bipolar NRZ');
subplot(4,1,2); stem(bnrz_filter005(47:mp*16+47));title('Bipolar NRZ fc = 0.05');
subplot(4,1,3); stem(bnrz_filter02(50:mp*16+50));title('Bipolar NRZ fc = 0.2');
subplot(4,1,4); stem(bnrz_filter04(50:mp*16+50));title('Bipolar NRZ fc = 0.4');
figure;
subplot(4,1,1); stem(xM(1:mp*16));title('Manchester');
subplot(4,1,2); stem(xM_filter005(47:mp*16+47));title('Manchester fc = 0.05');
subplot(4,1,3); stem(xM_filter02(50:mp*16+50));title('Manchester fc = 0.2');
subplot(4,1,4); stem(xM_filter04(50:mp*16+50));title('Manchester fc = 0.4');