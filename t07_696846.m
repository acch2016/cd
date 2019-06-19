%% ejercicio 1 ejercicio 2
%% 1
clc; clear all;
fs = 8000; B = 1000; Rb = 2000;
% Rb = fs/mp %bps
mp = fs/Rb; %4
beta = 0; Tp = 1/Rb; D=10;
type = 'rc';
Ts=1/fs;
energy =Tp;
[p t] = rcpulse(beta,D,Tp,Ts,type,energy);
% plot(t,p);
% Ep = Ts*sum(p.*p) % energía del pulso 0.01   Ts es el differential
% mp = Tp/Ts
% muestras_totales  = D*Tp %
% tiempo_total = D*Tp %s
wvtool(p)
%% 2
clear all; clc
fs = 8000; 
Ts = 1/fs; 
beta = 0.2; r = beta; B = 1000; Rs = 2*B/(1+r); 
Rb = Rs; 
Tp = 1/Rb;
Rb = 1600 % asi mp es de 5  % mp = fs/Rb
Tp = 1/Rb;
% mp = Tp/Ts; 
% Rb = fs/mp;
D = 10;
type = 'rc';
energy = Tp;
[p t] = rcpulse(beta,D,Tp,Ts,type,energy);
% plot(t,p);
wvtool(p)
%% 3
clear all; clc
fs = 4000; 
Ts = 1/fs; 
beta = 0.8;
D = 6;
Rb = 2000; Tp = 1/Rb;
type = 'rc';
energy = Tp;
[p t] = rcpulse(beta,D,Tp,Ts,type,energy);
% plot(t,p);
wvtool(p)

%% ejercicio 3
clc;clear all;
beta=0; fs=1000; Tp=1/100; D=10;
type = 'srrc';
Ts=1/fs;
energy =Tp;
[p t] = rcpulse(beta,D,Tp,Ts,type,energy);
fvtool(p);
% plot(t,p)
e = Ts*p*p'

bits = [1 0 1 1 0 0 1 1 1 1 1]
sym = bits*2-1
mp = Tp/Ts;
s = zeros(1,numel(bits)*mp);
s(1:mp:end) = sym;
xnrz = conv(p, s);
xnrz_Rx = conv(p, xnrz)*(1/mp);
plot(xnrz_Rx)
hold on;
stem([zeros(1,101) s])
%% señal Rx
nrz_rx = filter(flipr(p),1,nrz_tx)*(1/mp); %Match filter
start = filter_delay + 1;
% sampling_rx = 
%% bits lena
clc;clear all;
load lena512.mat
% whos('-file','lena512.mat')
% imshow(uint8(lena512))
lenarec=lena512(252:284,318:350);
% imshow(uint8(lenarec))
b=de2bi(lenarec,8);
% b=de2bi(lenarec,8,'left-msb'); % % % % % % % % % % % % % % version profesor
b=b';
bits=b(:);
%% pregunta 3

beta=0; fs=1000; Tp=1/100; D=10;
type = 'srrc';
Ts=1/fs;
energy =Tp;
[p t] = rcpulse(beta,D,Tp,Ts,type,energy);
fvtool(p);
% plot(t,p)
e = Ts*p*p'

sym = bits*2-1;
mp = Tp/Ts;
s = zeros(1,numel(bits)*mp);
s(1:mp:end) = sym;
xnrz = conv(p, s);
xnrz_Rx = conv(p, xnrz)*(1/mp);
plot(xnrz_Rx)
hold on;
stem([zeros(1,101) s])
%%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % reconstruir imagen
bR = reshape(bits_Rx, [8,1089]);
% isequal(bR,b)
bR = bR';
lena = reshape(bi2de(bR), size(lenarec));
% isequal(lenarec,lena)
figure; imshow(uint8(lena)); title([num2str(linecode), num2str(fc)]);
disp(['Errores: ', num2str( sum(xor(bits_Rx, bits')) ) ]);
disp('% % % % % EOF');
