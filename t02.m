clc; clear all;

% leer solo 10 segundos
% samples = [1, 10*44100];
% [y, Fs] = audioread('TheLovecats.wav',samples);
% audiowrite('TheLovecats441000.wav',y,Fs,'BitsPerSample',16);
% num2 = numel(y)

% leer el archivo ya cortado previamente
[y, Fs] = audioread('TheLovecats441000.wav');
info = audioinfo('TheLovecats441000.wav')
num = numel(y)
% y
miny = min(y)
maxy = max(y)

for index = 4:2:14
    % N puede tomar los siguientes valores
    % 2^14, 2^10, 2^8, 2^4
    % 16384, 1024, 256, 16
    N = 2^index;
    % paso de cuantizacion -> delta
    A = 1;
    delta = 2*A / N;
    % muestra cuantizada -> n
    xq = delta*(round(y/delta));
    if 4 == index
        audiowrite('C:\Users\ie696846\Documents\MATLAB\ie696846\TheLovecats16.wav',xq,Fs,'BitsPerSample',16);
        % audiowrite('/home/acc/Documents/MATLAB/ComunicacionesDigitales/TheLovecats4.wav', xq, Fs,'BitsPerSample',16);
        [y, Fs] = audioread('C:\Users\ie696846\Documents\MATLAB\ie696846\TheLovecats16.wav');
        % [y, Fs] = audioread('/home/acc/Documents/MATLAB/ComunicacionesDigitales/TheLovecats4.wav');
    end
    if 8 == index
        audiowrite('C:\Users\ie696846\Documents\MATLAB\ie696846\TheLovecats256.wav',xq,Fs,'BitsPerSample',16);
        % audiowrite('/home/acc/Documents/MATLAB/ComunicacionesDigitales/TheLovecats256.wav', xq, Fs,'BitsPerSample',16);
        [y, Fs] = audioread('C:\Users\ie696846\Documents\MATLAB\ie696846\TheLovecats256.wav');
        % [y, Fs] = audioread('/home/acc/Documents/MATLAB/ComunicacionesDigitales/TheLovecats256.wav');
    end
    if 10 == index
        audiowrite('C:\Users\ie696846\Documents\MATLAB\ie696846\TheLovecats1024.wav',xq,Fs,'BitsPerSample',16);
        % audiowrite('/home/acc/Documents/MATLAB/ComunicacionesDigitales/TheLovecats1024.wav', xq, Fs,'BitsPerSample',16);
        [y, Fs] = audioread('C:\Users\ie696846\Documents\MATLAB\ie696846\TheLovecats1024.wav');
        % [y, Fs] = audioread('/home/acc/Documents/MATLAB/ComunicacionesDigitales/TheLovecats1024.wav');
    end
    if 14 == index
        audiowrite('C:\Users\ie696846\Documents\MATLAB\ie696846\TheLovecats16384.wav',xq,Fs,'BitsPerSample',16);
        % audiowrite('/home/acc/Documents/MATLAB/ComunicacionesDigitales/TheLovecats16384.wav', xq, Fs,'BitsPerSample',16);
        [y, Fs] = audioread('C:\Users\ie696846\Documents\MATLAB\ie696846\TheLovecats16384.wav');
        % [y, Fs] = audioread('/home/acc/Documents/MATLAB/ComunicacionesDigitales/TheLovecats16384.wav');
    end
end

audiowrite('C:\Users\ie696846\Documents\MATLAB\ie696846\TheLovecats16.wav',xq,Fs,'BitsPerSample',16);
% audiowrite('/home/acc/Documents/MATLAB/ComunicacionesDigitales/TheLovecats4.wav', xq, Fs,'BitsPerSample',16);
[y, Fs] = audioread('C:\Users\ie696846\Documents\MATLAB\ie696846\TheLovecats16.wav');
% [y, Fs] = audioread('/home/acc/Documents/MATLAB/ComunicacionesDigitales/TheLovecats4.wav');
audio = audioplayer(y,Fs);
play(audio);

figure('Position',[100, 100, 400, 200]);
pspectrum(xq,Fs,'spectrogram','TimeResolution',0.2, 'OverlapPercent',99,'Leakage',0.85)




