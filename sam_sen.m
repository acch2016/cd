% Funcion que simula el muestreo de una se�al sinusoidal
% Formato:  [xa,xs]=sam_sen(Fo,Fs,Tw)
%           Fo es la frecuencia de operaci�n
%           Fs es la frecuencia de muestreo
%           Tw ventana de tiempo en milisegundos
%           N  es el n�mero de muestras de la secuencia

function [xa,xs] = sam_sen(Fo,Fs,Tw)
t=linspace(0,Tw*10^-3,10000);
xa=cos(2*pi*Fo*t);
% Se�al anal�gica compuesta de senoides pesadas
%xa=cos(2*pi*30*t) + cos(2*pi*150*t) + cos(2*pi*170*t)+ cos(2*pi*250*t) +cos(2*pi*330*t);
% Muestreo de la se�al 
N= round((Tw*10^-3)*Fs); % N�mero de muestras en la ventana de tiempo Tw
n=0:N;              
T=1/Fs;
nT= n*T;            % Instantes discretos tn=nT
wo=2*pi*Fo/Fs;      % Frecuencia angular digital normalizada
xs=cos(2*pi*Fo*nT);% Se�al muestrada; igual a xs=cos(wo*n);   
%xs=cos(2*pi*30*nT) + cos(2*pi*150*nT) + cos(2*pi*170*nT)+ cos(2*pi*250*nT) +cos(2*pi*330*nT);
figure(1);
plot(t,xa);         % Graficar la se�al anal�gica
grid on
hold on
stem(nT,xs,'filled','r');        % Graficar la se�al discretizada
title('Gr�fica para simular el muestreo');
xlabel('Tiempo')
ylabel('Magnitud')


% Omar L. rev. 22/02/2011