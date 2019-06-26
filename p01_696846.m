clc
clear

%% leer 10 segundos
% samples = [1, 10*44100];
% [y, Fs] = audioread('sippur.wav',samples);
% info = audioinfo('sippur.wav')
% audiowrite('/home/acc/Documents/MATLAB/ComunicacionesDigitales/sippur.wav', y, Fs,'BitsPerSample',16);
%% convertir el audio a mono
[y, Fs] = audioread('sippur.wav');
info = audioinfo('sippur.wav')
yizq = y(:,1);
yder = y(:,2);
y = (yizq+yder)/2; %convertir el audio a mono
subplot(2,1,1); pwelch(y,[],[],[],Fs,'power'); title('PSD of original signal');
% audio = audioplayer(y,Fs);
% play(audio);
%% filtrar la señal a 15kHz
m = [1 1 0 0];
o = 100;
B = 15e3;
fc = B/(Fs/2);
f = [0 fc fc 1];
LPF = fir2(o,f,m);
% wvtool(fil)
% fvtool(fil)
y = filter(LPF,1,y);
subplot(2,1,2); pwelch(y,[],[],[],Fs,'power'); title('PSD of filtered signal');
miny = min(y)
maxy = max(y)
% audio = audioplayer(y,Fs);
% play(audio);
% numel(y)
% numel(yf)
%% Normalizar la potencia de la señal filtrada a 1 watt
py = sum(y.^2)/numel(y)
y = y/sqrt(py);
% py = sum(y.^2)/numel(y)
% soundsc(y,Fs);
% play(audio);
%% Encontrar la potencia del ruido a la salida del filtro receptor, para N0 = 1/(15000*10^(0:0.3:3))
step = 0:0.3:3;
N0 = 1./(B*10.^step);
PN = B.*N0;
SNR = 1./PN;
SNRdB = 10*log10(SNR);
yorig = y;
for i = 1:numel(PN)
    Noise = sqrt(PN(i))*randn(1,numel(y));
%     minN = min(Noise)
%     maxN = max(Noise)
    y = yorig+Noise';
    y = filter(LPF,1,y);
    y = y./max(abs(y));
%     y = rescale(y,-1,1);
    audiowrite(['/home/acc/Documents/MATLAB/ComunicacionesDigitales/sippur', num2str(round(SNRdB(i))), '.wav'], y, Fs,'BitsPerSample',16);
end
% numel(Noise)
%%
fileID = fopen('/home/acc/Documents/MATLAB/ComunicacionesDigitales/sippur.wav')
file = fread(fileID)