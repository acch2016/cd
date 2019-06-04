clear all;
clc;
[y, Fs] = audioread('TheLovecats441000.wav');
% [y,Fs] = audioread('spring2s.wav');
% y es una matriz que tiene valores de las muestras las cuales de muestrear
% Fs = 44100
s = size(y);
s % 736556 filas
Y = fft(y);
dv = 44100/numel(y);
v = 0 : dv : 44100-dv;
figure;
plot(v,abs(Y));
%semilogx(v,abs(Y));
grid on;
xlim([0 44100]);% Set del límite del eje de las ordenadas
%text(1000,2500,' \leftarrow 1kHz','FontSize',18);
text(10000,2500,' \leftarrow 10kHz','FontSize',18);

fhi = [0 10000/44100 10000/44100 1];% en 1kHz no se percibe así que se lo puse a 10kHz
mhi = [0 0 1 1];
flo = [0 1000/44100 1000/44100 1];
mlo = [1 1 0 0];
for n=5:10;
    if n > 5 && n < 10;
        continue
    else
        bhi=fir2(n,fhi,mhi);
        %bhi
        blo=fir2(n,flo,mlo);
        %blo
        %Grafica de la respuesta en frecuencia del filtro
        figure;
        freqz(bhi,1)
        grid on;
        if n==5;
            title('HPF H(w) n=5','FontName','Courier','FontSize',15);
        else
            title('HPF H(w) n=10','FontName','Courier','FontSize',15);
        end
        figure;
        freqz(blo,1)
        grid on;
        if n==5;
            title('LPF H(w) n=5','FontName','Courier','FontSize',15);
        else
            title('LPF H(w) n=10','FontName','Courier','FontSize',15);
        end
    end
end

a=1;
yhi = filter(bhi,a,y); %out HPF
ylo = filter(blo,a,y); %out LPF
n = (0:length(y)-1)/Fs;

figure;
subplot(3,1,1)
plot(n,y)
grid on;
title('Original Signal','FontName','Courier','FontSize',15)
ys = ylim;

subplot(3,1,2)
plot(n,yhi)
grid on;
title('High Pass Filtered Signal','FontName','Courier','FontSize',15)
xlabel('n')
ylim(ys)

subplot(3,1,3)
plot(n,ylo)
grid on;
title('Low Pass Filtered Signal','FontName','Courier','FontSize',15)
xlabel('n')
ylim(ys)

audio=audioplayer(y,Fs);
%play(audio);
%pause;
audioHPF=audioplayer(yhi,Fs);
% play(audioHPF);
% pause;
audioLPF=audioplayer(ylo,Fs);
%play(audioLPF);
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                                   %TB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fp=800;
fr=1200;
% Cálculo de las frecuencias del filtro digital
wp = 2*pi*fp/Fs; % Frecuencia límite de banda de paso
wr = 2*pi*fr/Fs; % Frecuencia límite de banda de rechazo

Rp = 1; % dB % Rizo de banda de paso
Rr = 3; % dB % Atenuación de banda de rechazo

% Lo
fc = 900;% Frecuencia deseada de operación del filtro digital, en Hertz
%Wc_digital = (2*pi*fc)/Fs;% Frecuencia deseada de operación del filtro digital, en rad/muestra ( entre 0 y pi )
Ts = 1/Fs;% Periodo
wpana = (2/Ts)*tan(wp/2); % Frecuencia analógica en rad/seg
wrana = (2/Ts)*tan(wr/2);

[N,Wn] = buttord(wpana, wrana, Rp, Rr, 's');
%LPF
[Bs,As] = butter(N,Wn,'low','s'); % obtener coeficientes en s
w=logspace(2,5,500);
figure;freqs(Bs,As,w); % graficar coeficientes en s
title('Coeficientes en s LPF método TB','FontName','Courier','FontSize',15);
[Bz,Az] = bilinear(Bs,As,Fs);% Transformación de los coeficientes en z
figure;freqz(Bz,Az); % graficar coeficientes en z
title('Coeficientes en z LPF método TB','FontName','Courier','FontSize',15);
ytblo = filter(Bz,Az,y); % filtrar
figure;plot(v,abs(fft(ytblo))); % graficar señal filtrada
title('Respuesta en frecuencia LPF método TB','FontName','Courier','FontSize',15);
%HPF
[Bs,As] = butter(N,Wn,'high','s'); % obtener coeficientes en s
w=logspace(2,5,500);
figure;freqs(Bs,As,w); % graficar coeficientes en s
title('Coeficientes en s HPF método TB','FontName','Courier','FontSize',15);
[Bz,Az] = bilinear(Bs,As,Fs);% Transformación de los coeficientes en z
figure;freqz(Bz,Az); % graficar coeficientes en z
title('Coeficientes en z HPF método TB','FontName','Courier','FontSize',15);
ytblo = filter(Bz,Az,y); % filtrar
figure('Color','white');plot(v,abs(fft(ytblo))); % graficar señal filtrada
title('Respuesta en frecuencia HPF método TB','FontName','Courier','FontSize',15);

% pause;
% audio_y_TB_lo=audioplayer(ytblo,Fs);
% play(audio_y_TB_lo);