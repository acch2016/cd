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
% yizq = y(:,1);
% yder = y(:,2);
% y = (yizq+yder)/2; %convertir el audio a mono
% audiowrite('/home/acc/Documents/MATLAB/ComunicacionesDigitales/sippur.wav', y, Fs,'BitsPerSample',16);
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
py = sum(y.^2)/numel(y)
% soundsc(y,Fs);
% play(audio);
%% Encontrar la potencia del ruido a la salida del filtro receptor, para N0 = 1/(15000*10^(0:0.3:3))
step = 0:0.3:3;
N0 = 1./(B*10.^step);
PN = B.*N0;
SNR = 1./PN;
SNRdB = 10*log10(SNR);
%%
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
%% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % PARTE DIGITAL
% % % % % % % % % se puede correr el codigo sin ejecutar la parte analogica
% % % % % % % % %
clc;clear
Fs = 44100;
m = [1 1 0 0];
o = 100;
B = 15e3;
fc = B/(Fs/2);
f = [0 fc fc 1];
LPF = fir2(o,f,m);
%% Convert to binary version 1
% fileID = fopen('/home/acc/Documents/MATLAB/ComunicacionesDigitales/sippur.wav');
% file = fread(fileID);
% fbit = de2bi(file(:),8);
%% Convert to binary version 2
% wavdata = audioread('/home/acc/Documents/MATLAB/ComunicacionesDigitales/sippur.wav');
% b = de2bi( typecast(single(wavdata), 'uint16'), 16  );
% b=b';
% bits=b(:);
%% Convert to decimal (version 2 Recover test)
% bR = reshape(bits,size(b));
% % isequal(bR,b)
% bR = bR';
% wavdata2 = reshape( typecast( bi2de(bR), 'single' ), size(wavdata) );
%% Convert to binary version without typecast 
x = audioread('/home/acc/Documents/MATLAB/ComunicacionesDigitales/sippur.wav');
b = 16; %bits per sample
swing = (2^b-1)/2;
xq_int = round(x*swing+swing);
xq_bin = de2bi(xq_int,b,'left-msb');
xq_bin = xq_bin';
bits = xq_bin(:);
%% Design SRRC
beta = 0.35;
B = 15e3;
Rb = 2*B /(1+beta);
fs = 96000;
mp = ceil(fs/Rb);
Rb = fs/mp; % bps
D = 6;
Tp=1/Rb;
Ts=1/fs;
type = 'srrc';
energy = 1;
[p t] = rcpulse(beta,D,Tp,Ts,type,energy);
% wvtool(p)
% fvtool(p)
% plot(t,p)
e = Ts*p*p'
%% tren de pulsos pnrz
% sym = double(bits)*2-1; % symbols
sym = bits*2-1; % symbols
s = zeros(1,numel(bits)*mp);
s(1:mp:end) = sym; % tren de impulsos
pnrz = conv(p, s); % tren de pulsos % Shaping Pulse
% figure; plot(pnrz(1+80:mp*16+80)); title('polar NRZ');
%% Normalice el tren de pulsos para que tenga potencia unitaria
py = sum(pnrz.^2)/numel(pnrz) % 656.0370
pnrz = pnrz/sqrt(py);
py = sum(pnrz.^2)/numel(pnrz)
%% 
step = 0:0.3:3;
N0 = 1./(B*10.^step);
PN = B.*N0;
SNR = 1./PN;
SNRdB = 10*log10(SNR);

pnrz_o = pnrz;
for i = 1:numel(PN)
%     7 a. Genere un vector de ruido de la potencia correspondiente y de tamaño igual al tren de pulsos generado
    Noise = sqrt(PN(i))*randn(1,numel(pnrz_o));
    %     minN = min(Noise)
    %     maxN = max(Noise)
%     7 b. Al tren de pulsos generado añadale el ruido AWGN
    pnrz = pnrz_o+Noise;
%     8. Filtre el tren de pulsos más ruido con el LPF en el receptor
    pnrz = conv(LPF,pnrz);
%     9. Pase la salida del filtro LPF del receptor por el filtro acoplado
    pnrz = conv(p, pnrz); %match filter
%     realice el muestreo a la salida del filtro receptor
    y_s = pnrz(81:mp:end);
%     umbral de desicion e instante de observacion
%     numel(y_s)
    y_s = y_s(1:numel(bits));
    sym_Rx = sign(y_s);
    bits_Rx = (sym_Rx +1)/2; % bits_Rx es lo que se va a reconstruir
%     11. Reconstruccion
    Mbits = reshape( bits_Rx', size(xq_bin) );
    Mint = reshape( bi2de(Mbits', 'left-msb'), size(x) );
    Mint = Mint';    
%     minN = min(Mint)
%     maxN = max(Mint)
%     Mint = Mint./max(abs(Mint));
    wavdata2 = rescale(Mint,-1,1);
    audiowrite(['/home/acc/Documents/MATLAB/ComunicacionesDigitales/Dsippur', num2str(round(SNRdB(i))), '.wav'], wavdata2, Fs,'BitsPerSample',16);
end