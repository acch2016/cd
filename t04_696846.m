clear all; clc
mp = 9; % muestras por pulso
pbase = triang(mp); % pulso base triangular de mp muestras
t = 0:8;
b = [1 0 1 0 1 1 0 0];
stem(t,pbase)
ylim([0 1])
s = [1 0 0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 1];
s1 = b;
s1(s1==0) = -1;
s = zeros(1,(numel(b)-1)*mp+1);
s(1:mp:end) = s1;
x = conv(s,pbase);
% stem(x)
%% 2da sección
load lena512.mat
% whos('-file','lena512.mat')
% imshow(uint8(lena512))
lenarec=lena512(252:284,318:350);
% imshow(uint8(lenarec))
b=de2bi(lenarec,8);
% b=de2bi(lenarec,8,'left-msb'); % % % % % % % version profesor
b=b';
bits=b(:);
%% Realice y codifique el proceso inverso (bits a pixeles) para reconstruir la imagen.
bR = reshape(bits,[8,1089]);
% isequal(bR,b)
bR = bR';
lena = reshape(bi2de(bR),size(lenarec));
isequal(lenarec,lena)
imshow(uint8(lena))
% % % % % % % % % % % % % % % % % % % % % % version profesor
% Mbits = reshape(bits,8,33*33);
% vint = bi2de(Mbits','left-msb');
% Mint = vec2mat(vint,33,33);
% Mint(1,1)
% Mint = Mint';
% imshow(uint8(Mint))
%%  Unipolar NRZ % % % NOTA correr 2da seccion para obtener los bits
mp = 10; % samples per pulse
Fs = 96000;
Ts = 1/Fs;
% The bit rate is Rb= Rs= Fs / mp, because 1 bit= 1 %symbol and % every symbol has mpsamples per bit
Rs = Fs / mp;
% unipolar nrz
pbase = rectwin(mp); % Pulso completo;+
wvtool(pbase)
s = zeros(1,numel(bits)*mp);
s(1:mp:end) = bits;
subplot(3,1,1); stem(s(1:mp*16)); title('tren de impulsos');
unrz = conv(pbase,s);
subplot(3,1,3);pwelch(unrz,[],[],[],Fs,'power'); title('PSD of Unipolar NRZ');
subplot(3,1,2);plot(unrz(1:mp*16)); title('Uniolar NRZ');
%% Polar NRZ % % % NOTA correr 2da seccion para obtener los bits
mp = 10; % samples per pulse
Fs = 96000;
Ts = 1/Fs;
% The bit rate is Rb= Rs= Fs / mp, because 1 bit= 1 %symbol and % every symbol has mpsamples per bit
Rs = Fs / mp;
% polar nrz
pbase = rectwin(mp); % Pulso completo;
wvtool(pbase);
sym = (bits * 2)-1;
s = zeros(1,numel(bits)*mp);
s(1:mp:end) = sym;
subplot(3,1,1); stem(s(1:mp*16)); title('tren de impulsos');
pnrz = conv(pbase,s);
subplot(3,1,3);pwelch(pnrz,[],[],[],Fs,'power'); title('PSD of Polar NRZ');
subplot(3,1,2);plot(pnrz(1:mp*16)); title('Polar NRZ');
%% Bipolar NRZ % % % NOTA correr 2da seccion para obtener los bits
mp = 10; % samples per pulse
Fs = 96000;
Ts = 1/Fs;
% The bit rate is Rb= Rs= Fs / mp, because 1 bit= 1 %symbol and % every symbol has mpsamples per bit
Rs = Fs / mp;
% bipolar nrz
pbase = rectwin(mp); % Pulso completo;
wvtool(pbase)
%
am = mod(1:length(bits(bits == 1)), 2); % 0s y 1s del tamaño de los bits en 1 en la variable bits
am(am == 0) = -1; % alternate mark (vector de -1s y 1s)
bits(bits == 1) = am; % alternate mark inversion AMI (tres estados)
%
s = zeros(1,numel(bits)*mp);
s(1:mp:end) = bits;figure;
subplot(3,1,1);stem(s(1:mp*16)); title('tren de impulsos');
bnrz = conv(pbase,s);
subplot(3,1,3);pwelch(bnrz,[],[],[],Fs,'power'); title('PSD of Bipolar NRZ');
subplot(3,1,2);plot(bnrz(1:mp*16));title('Bipolar NRZ');
figure;powerbw(bnrz,Fs) % para obtener el bw de -3dB
%% Manchester % % % NOTA correr 2da seccion para obtener los bits
mp = 10; % samples per pulse
Fs = 96000;
Ts = 1/Fs;
% The bit rate is Rb= Rs= Fs / mp, because 1 bit= 1 %symbol and % every symbol has mpsamples per bit
Rs = Fs / mp
Rb = Rs
% Manchester
pM = [-ones(1,mp/2) ones(1,mp/2)]; %pulso base
wvtool(pM);
sym = (bits * 2)-1;
s = zeros(1,numel(bits)*mp);% tren de impulsos
s(1:mp:end) = sym; %tren
subplot(3,1,1); stem(s(1:mp*16)); title('tren de impulsos'); 
xM = conv(pM,s);
% subplot(3,1,3); pwelch(xM,200,100,200,Fs,'power'); title('PSD of Manchester');
subplot(3,1,3); pwelch(xM,[],[],[],Fs,'power'); title('PSD of Manchester');
subplot(3,1,2); plot(xM(1:mp*16));title('Manchester');
