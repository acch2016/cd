% mp = 9; % muestras por pulso
% pbase = triang(mp); % pulso base triangular de mp muestras
% t = 0:8;
% b = [1 0 1 0 1 1 0 0];
% stem(t,pbase)
% ylim([0 1])
% s = [1 0 0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 1];
% s1 = b;
% s1(s1==0) = -1;
% s = zeros(1,(numel(b)-1)*mp+1);
% s(1:mp:end) = s1;
% x = conv(s,pbase);
% stem(x)
load lena512.mat
% whos('-file','lena512.mat')
% imshow(uint8(lena512))
lenarec=lena512(252:284,318:350);
% imshow(uint8(lenarec))
b=de2bi(lenarec,8);
b=b';
bits=b(:);
%% Realice y codifique el proceso inverso (bits a pixeles) para reconstruir la imagen.
bR = reshape(bits,[8,1089]);
% isequal(bR,b)
bR = bR';
lena = reshape(bi2de(bR),size(lenarec));
isequal(lenarec,lena)
imshow(uint8(lena))
%%  Unipolar NRZ
% correr primera seccion para obtener los bits
mp = 10; % samples per pulse
Fs = 96000;
Ts = 1/Fs;
% The bit rate is Rb= Rs= Fs / mp, because 1 bit= 1 %symbol and % every symbol has mpsamples per bit
Rs = Fs / mp;
% unipolar nrz
pbase = rectwin(mp); % Pulso completo;
s = zeros(1,numel(bits)*mp);
s(1:mp:end) = bits;
stem(s(1:mp*16))
wvtool(pbase)
unrz = conv(pbase,s);
figure;pwelch(unrz,[],[],[],Fs,'power'); %title('PSD of Unipolar NRZ');
figure;plot(unrz(1:mp*16))
%% Polar NRZ
mp = 10; % samples per pulse
Fs = 96000;
Ts = 1/Fs;
% The bit rate is Rb= Rs= Fs / mp, because 1 bit= 1 %symbol and % every symbol has mpsamples per bit
Rs = Fs / mp;
% polar nrz
pbase = rectwin(mp); % Pulso completo;
sym = (bits * 2)-1
s = zeros(1,numel(bits)*mp);
s(1:mp:end) = sym;
stem(s(1:mp*16))
wvtool(pbase)
pnrz = conv(pbase,s);
figure;pwelch(pnrz,[],[],[],Fs,'power'); %title('PSD of Unipolar NRZ');
figure;plot(pnrz(1:mp*16))
%% Bipolar NRZ
mp = 10; % samples per pulse
Fs = 96000;
Ts = 1/Fs;
% The bit rate is Rb= Rs= Fs / mp, because 1 bit= 1 %symbol and % every symbol has mpsamples per bit
Rs = Fs / mp;
% bipolar nrz
pbase = rectwin(mp); % Pulso completo;
%
am = mod(1:length(bits(bits == 1)), 2); % 0s y 1s del tamaño de los bits en 1 en la variable bits
am(am == 0) = -1; % alternate mark (vector de -1s y 1s)
bits(bits == 1) = am; % alternate mark inversion AMI (tres estados)
%
s = zeros(1,numel(bits)*mp);
s(1:mp:end) = bits;
stem(s(1:mp*16))
wvtool(pbase)
pnrz = conv(pbase,s);
figure;pwelch(pnrz,[],[],[],Fs,'power'); %title('PSD of Unipolar NRZ');
figure;plot(pnrz(1:mp*16))
%% Manchester
mp = 10; % samples per pulse
Fs = 96000;
Ts = 1/Fs;
% The bit rate is Rb= Rs= Fs / mp, because 1 bit= 1 %symbol and % every symbol has mpsamples per bit
Rs = Fs / mp;
% Manchester
% pbase = rectwin(mp); % Pulso completo ???
pM = [-ones(1,mp/2) ones(1,mp/2)] %pulso base
sym = (bits * 2)-1
s = zeros(1,numel(bits)*mp);% tren de impulsos
% s(1:mp:end) = bits;
s(1:mp:end) = sym;%tren
stem(s(1:mp*16))
wvtool(pM)
xM = conv(pM,s);
figure;pwelch(xM,[],[],[],Fs,'power'); %title('PSD of Unipolar NRZ');
plot(xM(1:mp*16))