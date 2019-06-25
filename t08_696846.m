%%  Ejercicio 1 Diseñe un pulso RC p con el que vamos a transmitir
clear all;clc
% NOTA importante: beta > 0 para que realmente sea un coseno elevado
beta=0.25; fs=8000; D=10; Rb=1000; %bps
Tp=1/Rb;
Ts=1/fs;
mp = Tp/Ts;
type = 'srrc';
energy = 0.1;
[p t] = rcpulse(beta,D,Tp,Ts,type,energy);
% wvtool(p);
% fvtool(p);
% plot(t,p) %duracion total 0.01 %duracion parte central 0.002
e = Ts*p*p'
%%  Ejercicio 2 Diseñe cuatro filtros para usarse como candidatos a filtro receptor
disp('elegir pulso base: ');
option = input('1)acoplado / 2)pasabajas / 3)rectangular / 4)half-sine : ');
switch option
    case 1
        h = p;
    case 2
        fc = 0.1640625; % 0.1640625*fs/2 = 656.2500 Hz
        m = [1 1 0 0];
        o = numel(p); %81 % o sera? mp*D %80
        f = [0 fc fc 1];
        LPF = fir2(o,f,m);
        eLPF = Ts*sum(LPF.^2); % energía LPF
        LPF = LPF / sqrt(eLPF); % normalizar energia de LPF     
        h = sqrt(0.1)*LPF; % establecer energia de LPF a 0.1
        eLPF = Ts*sum(h.^2); % energía LPF
    case 3
        % muestras totales rc = D*mp %80
        % duracion total rc = D*Tp %0.01
        % rect = rectwin(2*mp);
        rect = ones(1,mp*2);
        eRect = Ts*sum(rect.^2); % energía rect
        rect = rect / sqrt(eRect); % normalizar energia de rect
        h = sqrt(0.1)*rect; % establecer energia de rect a 0.1
        efiltro = Ts*sum(h.^2); % energía rect
    case 4
        mb2 = round(mp*2);
        n = 0:mb2-1; w0=pi/(mb2);
        half_sin = sin(w0*n); % medio-seno 16 muestras
        ehalf_sin = Ts*sum(half_sin.^2); % energía half_sin
        half_sin = half_sin / sqrt(ehalf_sin); % normalizar energia de half_sin     
        h = sqrt(0.1)*half_sin; % establecer energia de half_sin a 0.1
        efiltro = Ts*sum(h.^2); % energía half_sin
end
% wvtool(filtro)
% fvtool(filtro)
%% Ejercicio 3. Encuentre la potencia del ruido a la salida de cada filtro
N = zeros(6,1);
for i = 1:6
    switch i
        case 1
            N02 = 0.1;
        case 2
            N02 = 0.08;
        case 3
            N02 = 0.06;
        case 4
            N02 = 0.04;
        case 5
            N02 = 0.025;
        case 6
            N02 = 0.01;
    end
    N(i) = N02*sum(h.*h)*Ts; % Noise Power % Observe que N es un escalar
end
%% ejercicio 4 
% Grafique la salida de cada uno de los filtros cuando la entrada es el pulso p
% Encuentre el valor maximo (con el comando max)
y = conv(p,h)*Ts; 
% figure('Position', [100, 500, 700, 200]); plot(y); 
% grid on;
% title(['Salida filtro receptor ', num2str(option)]);
% ylabel('Amplitud')
% xlabel('Muestras')
yMax = max(conv(p,h)*Ts)
%% ejercicio 5 SNR
SNR = zeros(6,1);
SNRdB = zeros(6,1);
% noise = zeros(6,1);
for i = 1:6
    SNR(i) = yMax^2/N(i);
    SNRdB(i) = 10*log10(SNR(i));
end
%% Ejercicio 6
counter = 0;
BER = zeros(6,1);
for i = 1:6
    noise = sqrt(N(i))*randn(1,1e6); % ruido AWGN
    for index = 1:numel(noise)
        
        if abs(noise(index)) > yMax
            counter = counter+1;
        end
    end
    BER(i) = (counter/2)/1e6;
    counter = 0;
    
end
% version eficiente
% (yMax-abs(noise)) < 0
% enseguida del vector resultante seguiria sumar los negativos
%% Ejercicio 7
% Grafique el BER para cada caso, SNR (eje horizontal), BER (en escala logarítmica)
% en el eje vertical
% figure('Position', [100, 500, 700, 200]);
% hold on;
semilogy(SNRdB, BER); %legend(['filtro ', num2str(option)])
grid on;
%% ejercicio 9 
load lena512.mat
% whos('-file','lena512.mat')
% imshow(uint8(lena512))
lenarec=lena512(252:284,318:350);
% imshow(uint8(lenarec))
b=de2bi(lenarec,8);
% b=de2bi(lenarec,8,'left-msb'); % % % % % % % % % % % % % % version profesor
b=b';
bits=b(:);
%% tren de pulsos pnrz
sym = bits*2-1;
s = zeros(1,numel(bits)*mp);
s(1:mp:end) = sym; % tren de impulsos
pnrz = conv(p, s); % tren de pulsos
%% c) Grafique la señal transmitida
time = 0 : Ts : Ts*numel(pnrz(1:16*mp+40))-Ts;
figure('Position', [100, 500, 700, 300]);
plot(time, pnrz(1:16*mp+40)/mp);
grid on;
title('Tx Signal')
ylabel('Amplitude')
xlabel('Time (s)')
hold on;
stem(time, [zeros(1,40) s(1:16*mp)])
%% Recibir
pnrz_Rx = conv(p, pnrz)*(1/mp); % receive with the Match Filter
time = 0 : Ts : Ts*numel(pnrz_Rx(1:16*mp+80))-Ts;
figure('Position', [100, 500, 700, 300]);
plot(time, pnrz_Rx(1:16*mp+80)/(mp+80))
title('Rx Signal');
ylabel('Amplitude')
xlabel('Time (s)')
hold on;
stem(time, [zeros(1,80) s(1:16*mp)])
%% Muestree la señal recibida en esos instantes
% umbral de desicion e instante de observacion
y_s = pnrz_Rx(81:mp:end);
numel(y_s)
y_s = y_s(1:8712);
sym_Rx = sign(y_s);
bits_Rx = (sym_Rx +1)/2; % bits_Rx es lo que se va a reconstruir
%% scatterplot(sampled_signal)
scatterplot(y_s)
%% d) Obtenga el espectro del tren de pulsos utilizando un estimador espectral de potencia.
figure('Position', [100, 500, 700, 200]);
pwelch(bits_Rx,[],[],[],fs,'power'); title('PSD of Rx pnrz ');
%% e) Obtenga el diagrama de ojo de la señal transmitida para 3 UI (unit interval)
eyediagram(pnrz,3*mp)
%% Ejercicio 10 filtros LPF
m = [1 1 0 0];
o = 100;
fc = input('ingrese una fc (".2" o ".4"): ');
f = [0 fc fc 1];
fil = fir2(o,f,m);
% wvtool(fil)
% fvtool(fil)
rx = conv(fil,pnrz); % señal recibida
% y=conv(p,h)*Ts
%% b) diagramas de ojo
eyediagram(rx,3*mp)
%% c) pase la salida de los filtros anteriores por el MATCH FILTER & eyediagram
pnrz_fil_RX = conv(p, rx)*(1/mp); % receive with the Match Filter
eyediagram(pnrz_fil_RX,3*mp)
%% d) umbral de desicion e instante de observacion 
y_s2 = pnrz_Rx(81:mp:end);
numel(y_s2)
y_s2 = y_s2(1:8712);
sym_Rx2 = sign(y_s2);
bits_Rx2 = (sym_Rx2 +1)/2; % bits_Rx es lo que se va a reconstruir
%% recuperar imagen % % % % % % % % % % % % % % % % % % % % % % % % % % % % % reconstruir imagen
bR = reshape(bits_Rx2, [8,1089]);
% isequal(bR,b)
bR = bR';
lena = reshape(bi2de(bR), size(lenarec));
% isequal(lenarec,lena)
figure; imshow(uint8(lena)); %title([num2str(linecode), num2str(fc)]);
disp(['Errores: ', num2str( sum(xor(bits_Rx2, bits')) ) ]);
% disp('% % % % % EOF');
