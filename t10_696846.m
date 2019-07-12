%% 
clear
clc
%%
randi([0, 1], [1, 60e3]);
%% Diseño del pulso base RC
beta = input('ingrese una beta (0/.2/.5/.9): ');
Rb = 11025;
B = Rb*(1+beta);
fs = 22050;
mp = ceil(fs/Rb);
Rb = fs/mp; % bps
D = 12;
Tp=1/Rb;
Ts=1/fs;
type = 'rc';
energy = Tp;
[pe t] = rcpulse(beta,D*2,Tp/2,Ts,type,energy);
% e = Ts*pe*pe'
deltas = [1 zeros(1, mp/2) -1]; % stem(deltas)
pman = conv(deltas, pe); % stem(pman)
% wvtool(pman)
pbman = pman(4:64);
wvtool(pbman)