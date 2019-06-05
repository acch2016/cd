clc; clear all;

% leer solo 10 segundos
% samples = [1, 10*44100];
% [y, Fs] = audioread('TheLovecats.wav',samples);
% audiowrite('TheLovecats441000.wav',y,Fs,'BitsPerSample',16);
% num2 = numel(y)

% leer el archivo ya cortado previamente
[y, Fs] = audioread('TheLovecats441000.wav');
info = audioinfo('TheLovecats441000.wav')
muestras = numel(y) % comprobando 44100*10 muestras
% comprobando la amplitud pico a pico
miny = min(y)
maxy = max(y)
 
N = 16; % N puede tomar los valores 2^14=16384, 2^10=1024, 2^8=256, 2^4=16
delta = 2 / N; % paso de cuantizacion
xq = delta*(round(y/delta)); % muestra cuantizada
% audiowrite('C:\Users\ie696846\Documents\MATLAB\ie696846\TheLovecats16.wav',xq,Fs,'BitsPerSample',16);
% % audiowrite('/home/acc/Documents/MATLAB/ComunicacionesDigitales/TheLovecats16.wav', xq, Fs,'BitsPerSample',16);

% [y, Fs] = audioread('C:\Users\ie696846\Documents\MATLAB\ie696846\TheLovecats16.wav');
% [y, Fs] = audioread('/home/acc/Documents/MATLAB/ComunicacionesDigitales/TheLovecats16.wav');
% audio = audioplayer(y,Fs);
% play(audio);

% figure('Position',[100, 100, 400, 200]);
% pspectrum(xq,Fs,'spectrogram','TimeResolution',0.2, 'OverlapPercent',99,'Leakage',0.85)

Fm = 10e3;
sampleFreq = Fs
f = [0 Fm/Fs Fm/Fs 1];
m = [1 1 0 0];
o = 99;
b = fir2(o,f,m);
[h, w] = freqz(b); % freqz calcula la respuesta en frecuencia de b
figure;
subplot(2,1,1)
plot(f,m,w/pi,abs(h));
title('Filtro 10kHz')
ylabel('Mag espectro')
xlabel('Frecuencia normalizada (Hz)')
lgs = {'respuesta deseada','respuesta obtenida'};
legend(lgs)
subplot(2,1,2)
semilogy(f,m,w/pi,abs(h));
title('Filtro 10kHz')
ylabel('Mag espectro (dB)')
xlabel('Frecuencia normalizada (Hz)')
lgs = {'respuesta deseada','respuesta obtenida'};
legend(lgs)
x = conv(y,b); % filtrado de y con filtro b

Y = fft(y);
dv = 44100/numel(y);
v = 0 : dv : 44100-dv;
figure;
subplot(2,1,1)
% plot(v,abs(Y));
semilogx(v,abs(Y));
title('Espectro señal original')
ylabel('Mag espectro')
xlabel('Frecuencia (log)')

X = fft(x);
dv = 44100/numel(x);
v = 0 : dv : 44100-dv;
subplot(2,1,2)
% plot(v,abs(Y));
semilogx(v,abs(X));
title('Espectro señal filtrada')
ylabel('Mag espectro')
xlabel('Frecuencia (log)')

B = rescale(x,-1,1);
% audiowrite('C:\Users\ie696846\Documents\MATLAB\ie696846\TheLovecats_f10k.wav',x,Fs,'BitsPerSample',16);
% audiowrite('/home/acc/Documents/MATLAB/ComunicacionesDigitales/TheLovecats_f10k.wav', B, Fs,'BitsPerSample',16);

figure('Position',[100, 100, 400, 200]);
pspectrum(B,Fs,'spectrogram','TimeResolution',0.2, 'OverlapPercent',99,'Leakage',0.85)
