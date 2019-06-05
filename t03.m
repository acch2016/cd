Fs = 100; Ts = 1/Fs; t = 0:Ts:1; 
x = zeros(1,numel(t)); 
x(find( (t>0.4) & (t<0.6) )) = 1; 
plot(t,x)
title('Pulso cuadrado, duración = 0.2s');
% xlim([0 100])
% wvtool(x)
figure;
pspectrum(x)