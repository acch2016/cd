clc; clear all;
% leemos el original
[y, Fs] = audioread('TheLovecats441000.wav');
info = audioinfo('TheLovecats441000.wav')
num = numel(y)
% ahora leemos solo 10 muestras
% samples = [1,10*Fs];
% [y, Fs] = audioread('TheLovecats.wav',samples);
% audiowrite('TheLovecats441000.wav',y,Fs,'BitsPerSample',16);
% num2 = numel(y)
y
miny = min(y)
maxy = max(y)
% audio = audioplayer(y,Fs);
% play(audio);
% pause

% 14, 10, 8 y 4

% audiowrite('TheLovecats8.wav',y,Fs,'BitsPerSample',8);
% [y, Fs] = audioread('TheLovecats8.wav');
% audio = audioplayer(y,Fs);
% play(audio);
% info = audioinfo('TheLovecats8.wav')
% b = 16;
version = 2;
if (1 == version)
    % % % version sin pasar a bits  % Cuantizar a entero
    b = 4;
    swing = (2^b-1)/2 %14->8191.5 %10->512 %8->128 %4->8
    xq_i    nt =(y*swing);
    max(xq_int)
    min(xq_int)
    class(xq_int)
    B = rescale(xq_int,-1,1)
end

% %  https://la.mathworks.com/matlabcentral/answers/48153-encoding-the-signal-with-14-bits
% D = dec2bin( (xq_int), 14)
% D = dec2bin(yswing);
% D = D(:,1:14)

b = 4; % Cantidad de bits 
% % Cuantizar a entero y expresión en binario
swing = (2^b-1)/2; 
xq_int = round(y*swing+swing);
max_xq_int = max(xq_int)
min_xq_int = min(xq_int)
xq_bin = de2bi(xq_int,14,'left-msb')
decimal = bi2de(xq_bin,'left-msb') % Convert binary vectors to decimal numbers
% tipo = class(decimal)
% min_decimal = min(decimal)
% max_decimal = max(decimal)
% decimal-swing

% yswing = y*swing
% min(yswing)
% max(yswing)
% xq_int = round(y*swing);
% xq_int
% min(xq_int)
% max(xq_int)
% F1 audiowrite. xq_int es un double y si el valor esta fuera de -1 y +1 => clipping
audiowrite('C:\Users\ie696846\Documents\MATLAB\ie696846\TheLovecats4.wav',B,Fs,'BitsPerSample',16);
[y, Fs] = audioread('C:\Users\ie696846\Documents\MATLAB\ie696846\TheLovecats4.wav');
% numel_xq_int = numel(y)
audio = audioplayer(y,Fs);
play(audio);
% soundsc(y,Fs,16) % Normalizado por hardware entre 1 y -1 volts
% info = audioinfo('TheLovecats10.wav')


