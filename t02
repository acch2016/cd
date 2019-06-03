clc; clear all;
% leemos el original
[y, Fs] = audioread('TheLovecats.wav');
info = audioinfo('TheLovecats.wav')
num = numel(y)
% ahora leemos solo 10 muestras
samples = [1,10*Fs];
[y, Fs] = audioread('TheLovecats.wav',samples);
num2 = numel(y)
% y
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

% Cuantizar a entero
b = 4;
swing = (2^b-1)/2 %14->8191.5 %10->512 %8->128 %4->8
xq_int = round(y*swing+swing);
min(xq_int)
max(xq_int)
class(xq_int)
B = rescale(uint8(xq_int),-1,1)
% %  https://la.mathworks.com/matlabcentral/answers/48153-encoding-the-signal-with-14-bits
% D = dec2bin( (xq_int), 14)
% D = dec2bin(yswing);
% D = D(:,1:14)

xq_bin = de2bi(xq_int,14,'left-msb')
decimal = bi2de(xq_bin,'left-msb') % Convert binary vectors to decimal numbers
tipo = class(decimal)
min_decimal = min(decimal)
max_decimal = max(decimal)
decimal-swing

yswing = y*swing
min(yswing)
max(yswing)
xq_int = round(y*swing);
xq_int
min(xq_int)
max(xq_int)
%% F1 audiowrite. xq_int es un double y si el valor esta fuera de -1 y +1 => clipping
audiowrite('TheLovecats10.wav',xq_int,Fs,'BitsPerSample',16);
[y, Fs] = audioread('TheLovecats10.wav');
% numel_xq_int = numel(y)
% audio = audioplayer(y,Fs);
% soundsc(y,Fs,16) % Normalizado por hardware entre 1 y -1 volts
% info = audioinfo('TheLovecats10.wav')


