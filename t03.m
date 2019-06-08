%% Ejercicio 1 pregutas 1,2,3
Fs = 100; Ts = 1/Fs; t = 0:Ts:1; 
x = zeros(1,numel(t)); 
x(find( (t>0.4) & (t<0.6) )) = 1; 
figure;plot(t,x)
title('Pulso cuadrado, duracion = 0.2s');
wvtool(x)
% figure;pspectrum(x)
% figure;pwelch(x,101,30,500,'one-side','power',Fs)
% figure;pwelch(x,101,30,500,'centered','power',Fs)
%% Ejercicio 1 preguta 4
Fs = 100; Ts = 1/Fs; t = 0:Ts:1; 
x = zeros(1,numel(t));
x(find( t<0.2 )) = 1;
plot(t,x)
title('Pulso cuadrado, duracion = 0.2s');
wvtool(x)
%% Ejercicio 2
clc;clear all
format compact
Fs = 100; Ts = 1/Fs; t = 0:Ts:1; 
x = zeros(1,numel(t)); 
x(find( (t>0.2) & (t<0.4) )) = 1; 
x(find( (t>0.7) & (t<0.9) )) = 1; 
figure;plot(t,x)
title('2 pulsos cuadrados, duracion = 0.2s');
wvtool(x)
%% Ejercicio 3 
clc;clear all
format compact
Fs = 100; Ts = 1/Fs; t = 0:Ts:1; 
x = zeros(1,numel(t)); 
x(find( (t>0.2) & (t<0.4) )) = 1; 
x(find( (t>0.7) & (t<0.9) )) = -1; 
figure;plot(t,x)
title('2 pulsos cuadrados bipolares, duracion = 0.2s');
wvtool(x)
%% Ejercicio 4 pregunta 1
clc;clear all
format compact
Fs = 100; Ts = 1/Fs; t = 0:Ts:1; 
x = zeros(1,numel(t)); 
x(find( (t>0.45) & (t<0.55) )) = 1; 
figure;plot(t,x)
title('Pulso cuadrado, duracion = 0.1s');
wvtool(x)