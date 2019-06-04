% Energy Signal
Fs=10;      % Sampling Frequency
Ts=1/Fs;          % Sampling Time 
t = [0:.1:2*pi];  % Times at which to sample the sine function 
Wo= 1;            % 1 rad 
sigE = sin(Wo*t); % Original signal, a sine wave 
Energy = Ts*sum(sigE.^2); % Otra forma (sigE*sigE')*Ts; 
Rxx= Ts*xcorr(sigE,sigE); % Autocorrelation of Energy Signal % Ts*conv(sigE,fliplr(sigE)) 
plot(Rxx)
%% Power Signal 
sigP= randn(1,1e6); 
Power = (sigP*sigP')/numel(sigP); 
Rxx = xcorr(sigP,sigP); 
plot(Rxx)

plot(sigE)
Energy
freqz(Rxx)
pwelch(sigP,500,300,500,'one-side','power')
%%
%Grafique el espectro de un pulso cuadrado 
Fs = 100; Ts = 1/Fs; t = 0:Ts:1; x = zeros(1,length(t)); x(find( (t>=0.4) & (t<0.6) )) = 1; 
title('Pulso cuadrado, duración = 0.2s'); % Analice el Ancho de Banda (según criterios)

plot(x)
Energy = Ts*sum(x.^2) %0.2 Joules
figure;wvtool(x)
0.1*Fs/2
%%
Fs = 1000;   t = 0:1/Fs:0.2; x = 3*cos(2*pi*t*200); 
plot(t,x)
% Using a Welch spectral estimator.
pwelch(x,[],[],[],'one-side','power',Fs)

Fs = 1000;   t = 0:1/Fs:10; x = 3*cos(2*pi*t*200);
figure; pwelch(x,[],[],[],'one-side','power',Fs)