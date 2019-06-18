
clc; clear all;beta=0; fs=1000; Tp=1/100; D=10;
type = 'rc';
Ts=1/fs;
energy =Tp;
[p t] = rcpulse(beta,D,Tp,Ts,type,energy);
plot(t,p);
Ep = Ts*sum(p.*p) % energía del pulso 0.01   Ts es el differential
mp = Tp/Ts
Rb = fs/mp %bps
% Rs = 100 bauds
muestras_totales  = D*Tp %100
tiempo_total = D*Tp %0.1 segundos
wvtool(p)
%% no salio
clc; clear all;beta=0; fs=10; Tp=1; D=12;
type = 'rc';
Ts=1/fs;
energy =Tp;
[p t] = rcpulse(beta,D,Tp,Ts,type,energy);
plot(t,p);
Ep = Ts*sum(p.*p) % energía del pulso 0.01   Ts es el differential
mp = Tp/Ts
Rb = fs/mp %bps
% Rs = 100 bauds
muestras_totales  =D*Tp %100
tiempo_total = D*Tp %0.1 segundos
wvtool(p)
%%
beta=0; Fs=10; mp=10; D=12;
Ts=1/Fs; fN=Fs/2;
Rs=Fs/mp;
Rb=Rs;
Tp=1/Rb; % In this case Rb=1 bps
type='rc'; energy=Tp; % Initial energy
[p t] = rcpulse(beta,D,Tp,Ts,type,energy); %Pulse generation
e=Ts*sum(p.^2); % Energy computation of the discrete pulse
% e=Ts*max(xcorr(p)); % e= (Ts)*p*p' % e = trapz(t,p.*p);
p=p./sqrt(e);
% Pulse with Energy = 1
wvtool(p)
%%
%Consideraciones del pulso formador
beta=0.35; fs=1000; Tp=1/100; D=10; energy=1;
type='rc';Ts=1/fs;
[p t] = rcpulse(beta,D,Tp,Ts,type,energy);
e=Ts*max(xcorr(p)); % otra forma: e = trapz(t,r.*r)
p=p./sqrt(e); % Normalizar la energia=1
%%

clear all 
clc
% diseño del filtro
m = [1 1 0 0];
o = 100;
fc = input('ingrese una fc (.046/.146/.246/.8): ');
f = [0 fc fc 1];
fil = fir2(o,f,m);
% fvtool(fil)
load lena512.mat
% whos('-file','lena512.mat')
% imshow(uint8(lena512))
lenarec=lena512(252:284,318:350);
% imshow(uint8(lenarec))
b=de2bi(lenarec,8);
% b=de2bi(lenarec,8,'left-msb'); % % % % % % % % % % % % % % version profesor
b=b';
bits = b(:);
bits = [1 0 1 1 0 0 1 1 1 1 1];

mp = 10; % samples per pulse
Fs = 1000;
Ts = 1/Fs;
% The bit rate is Rb= Rs= Fs / mp, because 1 bit= 1 %symbol and % every symbol has mpsamples per bit
Rs = Fs / mp;

p = input('elegir pulso base (1 cuadrado / 2 medio-seno RZ / 3 sinc): ');
switch p
    case 1
        pbase = rectwin(mp); % Pulso completo UNRZ PNRZ BNRZ;
        pm = [ones(1,mp/2) -ones(1,mp/2)]; %pulso base Manchester opcion A
    case 2
        mb2 = round(mp/2);
        n = 0:mb2-1; w0=pi/(mb2);
        half_sin = sin(w0*n); % Probar con medio-seno RZ como en Zigbee
        pm = [half_sin -half_sin]; pbase = pm;
    case 3
        beta=0; fs=1000; Tp=1/100; D=10;
        type = 'rc';
        Ts=1/fs;
        energy =Tp;
        [p t] = rcpulse(beta,D,Tp,Ts,type,energy);
end

if linecode == 'pnrz'
    
    wvtool(pbase);
    sym = (bits * 2)-1;
    s = zeros(1,numel(bits)*mp);
    s(1:mp:end) = sym;
    subplot(3,1,1); stem(s(1:mp*16)); title('tren de impulsos');
    pnrz = conv(pbase,s);
    subplot(3,1,2); plot(pnrz(1:mp*16)); title('Polar NRZ');
    subplot(3,1,3); figure;pwelch(pnrz,[],[],[],Fs,'power'); title('PSD of Polar NRZ');    
    
    figure;
    subplot(2,2,1); stem(pnrz(1:mp*16)); title('Polar NRZ');
    subplot(2,2,2); pwelch(pnrz,[],[],[],Fs,'power'); title('PSD of Polar NRZ');
    rx = conv(fil,pnrz); % señal recibida
    subplot(2,2,3); stem(rx(1:mp*16)); title(['Polar NRZ fc = ', num2str(fc)]);
    subplot(2,2,4); pwelch(rx,[],[],[],Fs,'power'); title(['PSD of Polar NRZ fc = ', num2str(fc)]);
    
    % % % % % % % % % % % % % % % instantes de observación & umbral de desicion
    y_s = rx(55:mp:end);
%     numel(y_s) % 8718
    y_s = y_s(1:8712);
    % scatterplot(y_s_pnrz_f1)
    sym_Rx = sign(y_s);
    bits_Rx = (sym_Rx +1)/2; % bits_Rx es lo que se va a reconstruir
    
end



% % % % % % % % % % % % % % % % % % % % % % % % % % % % % reconstruir imagen
bR = reshape(bits_Rx, [8,1089]);
% isequal(bR,b)
bR = bR';
lena = reshape(bi2de(bR), size(lenarec));
% isequal(lenarec,lena)
figure; imshow(uint8(lena)); title([num2str(linecode), num2str(fc)]);
disp(['Errores: ', num2str( sum(xor(bits_Rx, bits')) ) ]);
disp('% % % % % EOF');


