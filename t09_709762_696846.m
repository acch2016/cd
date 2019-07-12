%% Diseño del pulso base SRRC original
clear
beta = 0.5;
B = 7.2e3;
Rb = 2*B /(1+beta);
fs = 96000;
mp = fs/Rb;
% Rb = fs/mp; % bps
D = 6;
Tp=1/Rb;
Ts=1/fs;
type = 'srrc';
energy = Tp;
[p t] = rcpulse(beta,D,Tp,Ts,type,energy);
wvtool(p)
%% Diseño del pulso base SRRC original 2
clear
beta = 0.25;
% B = 7.2e3;
% Rb = 2*B /(1+beta);
fs = 96000;
mp = 6;
Rb = fs/mp; % bps
B = Rb*(1+beta)/2;
D = 6;
Tp=1/Rb;
Ts=1/fs;
type = 'srrc';
energy = Tp;
[p t] = rcpulse(beta,D,Tp,Ts,type,energy);
wvtool(p)
%% Diseño del pulso base SRRC
clear
beta = 0.25; %.2 .35
% B = 7.2e3;
% Rb = 2*B /(1+beta);
fs = 96000;
mp = 6;
Rb = fs/mp; % bps
B = Rb*(1+beta)/2;
% mp = fs/Rb;
% Rb = fs/mp; % bps
D = 6;
Tp=1/Rb;
Ts=1/fs;
type = 'srrc';
energy = Tp;
[pe t] = rcpulse(beta,D*2,Tp/2,Ts,type,energy);
% e = Ts*pe*pe'
deltas = [1 zeros(1, mp/2) -1]; % stem(deltas)
pman = conv(deltas, pe); % stem(pman)
% wvtool(pman)
pbman = pman(3:39);
% wvtool(pbman)
%% LENA 
load lena512.mat
whos('-file','lena512.mat')
% imshow(uint8(lena512))
% lenarec=lena512(252:284,318:350);
% lenarec=lena512(242:284,308:350);
% lenarec=lena512(242:294,308:360);
lenarec=lena512(1:512,1:512);
% imshow(uint8(lenarec))
b=de2bi(lenarec,8);
% b=de2bi(lenarec,8,'left-msb'); % % % % % % % % % % % % % version profesor
b=b';
bits=b(:);
%% [ Preambulo Header Payload ]
pr = [1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 1 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 1];
altoancho =  de2bi(size(lenarec), 16, 'left-msb');
altoancho=altoancho';
header = altoancho(:);
% header(1,9:16) = altoancho(2,:);
% agregamos longitud del header
lengthheader = de2bi(length(header), 8, 'left-msb');
bits = [pr'; lengthheader'; header; bits];
%% tren de pulsos Manchester
sym = bits*2-1; % symbols
s = zeros(1,numel(bits)*mp); % zero vector with length of data
s(1:mp:end) = sym; % tren de impulsos
m = conv(pbman, s); %tren de pulsos
% figure('Position',[100, 100, 700, 250]); plot(m(1:30*mp));
% title('Transmitted signal')
% ylabel('Amplitude (V)')
% xlabel('Time (s)')
%% Normalize Power
pm = sum(m.^2)/numel(m) 
m = m/sqrt(pm);
pm = sum(m.^2)/numel(m)
%% PSD
figure('Position',[100, 100, 700, 250]); pwelch(m, 500,300,500,'one-side', 'power', fs);
%% Tx
soundsc([zeros(1, fs/2) m zeros(1, fs/2)],fs);