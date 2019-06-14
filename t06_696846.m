clear all; clc;
% diseño del filtro
m = [1 1 0 0];
o = 100;
fprintf('para facilitar la revision se desplegaran instrucciones \n para elegir el linecode y el filtro \n\r');
disp('filtro 1 fc = 0.046')
disp('filtro 2 fc = 0.146')
disp('filtro 3 fc = 0.246')
disp('filtro 4 fc = 0.546')
disp('filtro 5 fc = 0.8')
fc = input('seleccione numero de filtro')
switch fc
    case 1
        f = [0 0.046 0.046 1];
        f1 = fir2(o,f,m);
    case 2
        
end

f = [0 0.046 0.046 1];

% fvtool(f1)
f = [0 0.146 0.146 1];
f2 = fir2(o,f,m);
% fvtool(f2)
f = [0 0.246 0.246 1];
f3 = fir2(o,f,m);
% fvtool(f3)
f = [0 0.546 0.546 1];
f4 = fir2(o,f,m);
% fvtool(f4)
f = [0 0.8 0.8 1];
f5 = fir2(o,f,m);
% fvtool(f5)
%%
load lena512.mat
% whos('-file','lena512.mat')
% imshow(uint8(lena512))
lenarec=lena512(252:284,318:350);
% imshow(uint8(lenarec))
b=de2bi(lenarec,8);
% b=de2bi(lenarec,8,'left-msb'); % % % % % % % version profesor
b=b';
bits=b(:);
%% Unipolar NRZ
mp = 10; % samples per pulse
Fs = 96000;
Ts = 1/Fs;
% The bit rate is Rb= Rs= Fs / mp, because 1 bit= 1 %symbol and % every symbol has mpsamples per bit
Rs = Fs / mp;
pbase = rectwin(mp); % Pulso completo;
% wvtool(pbase)
s = zeros(1,numel(bits)*mp);
s(1:mp:end) = bits;
% subplot(3,1,1); stem(s(1:mp*16)); title('tren de impulsos');
unrz = conv(pbase,s);
% subplot(3,1,3);
% figure;pwelch(unrz,[],[],[],Fs,'power'); title('PSD of Unipolar NRZ');
%% Polar NRZ % 
mp = 10; % samples per pulse
Fs = 96000;
Ts = 1/Fs;
% The bit rate is Rb= Rs= Fs / mp, because 1 bit= 1 %symbol and % every symbol has mpsamples per bit
Rs = Fs / mp;
pbase = rectwin(mp); % Pulso completo;
% wvtool(pbase);
sym = (bits * 2)-1;
s = zeros(1,numel(bits)*mp);
s(1:mp:end) = sym;
% subplot(3,1,1); stem(s(1:mp*16)); title('tren de impulsos');
pnrz = conv(pbase,s);
% subplot(3,1,3);
% figure;pwelch(pnrz,[],[],[],Fs,'power'); title('PSD of Polar NRZ');
%% Bipolar NRZ % 
mp = 10; % samples per pulse
Fs = 96000;
Ts = 1/Fs;
% The bit rate is Rb= Rs= Fs / mp, because 1 bit= 1 %symbol and % every symbol has mpsamples per bit
Rs = Fs / mp;
pbase = rectwin(mp); % Pulso completo;
% wvtool(pbase)
am = mod(1:length(bits(bits == 1)), 2); % 0s y 1s del tamaño de los bits en 1 en la variable bits
am(am == 0) = -1; % alternate mark (vector de -1s y 1s)
bits1 = bits; % para no modificar la variable original que contine los bits
bits1(bits == 1) = am; % alternate mark inversion AMI (tres estados)
s = zeros(1,numel(bits1)*mp);
s(1:mp:end) = bits1;
% subplot(3,1,1);stem(s(1:mp*16)); title('tren de impulsos');
bnrz = conv(pbase,s);
% subplot(3,1,3);
% figure;pwelch(bnrz,[],[],[],Fs,'power'); title('PSD of Bipolar NRZ');
% subplot(3,1,2);plot(bnrz(1:mp*16));title('Bipolar NRZ');
% figure;powerbw(bnrz,Fs) % para obtener el bw de -3dB
%% Manchester % %
mp = 10; % samples per pulse
Fs = 96000;
Ts = 1/Fs;
% The bit rate is Rb= Rs= Fs / mp, because 1 bit= 1 %symbol and % every symbol has mpsamples per bit
Rs = Fs / mp;
pbase = rectwin(mp); % Pulso completo;
pM = [-ones(1,mp/2) ones(1,mp/2)]; %pulso base
wvtool(pM);
sym = (bits * 2)-1;
s = zeros(1,numel(bits)*mp);% tren de impulsos
s(1:mp:end) = sym; %tren
% subplot(3,1,1); stem(s(1:mp*16)); title('tren de impulsos'); 
xM = conv(pM,s);
% subplot(3,1,3); pwelch(xM,200,100,200,Fs,'power'); title('PSD of Manchester');
% subplot(3,1,3); 
% figure;pwelch(xM,[],[],[],Fs,'power'); title('PSD of Manchester');
% subplot(3,1,2); plot(xM(1:mp*16));title('Manchester');
%% filter 0.546 ejemplo en clase
unrz_filter0546 = conv(f4,unrz);
% pwelch(unrz_filter0546,[],[],[],Fs,'power'); title('PSD Unipolar NRZ fc = 0.546');
subplot(2,1,1); stem(unrz(1:mp*16)); title('Unipolar NRZ');
subplot(2,1,2); stem(unrz_filter0546(1:mp*16)); title('Unipolar NRZ fc = 0.546');

y_s = unrz_filter0546(55:mp:end);
numel(y_s) % 8718
y_s = y_s(1:8712);
% scatterplot(y_s)
sym_Rx = sign(y_s-0.5);
bits_Rx = (sym_Rx+1)/2;
%% Filtre la señal Unipolar NRZ con los tres filtros. La salida del filtro representa la señal recibida
figure;
subplot(4,2,1); stem(unrz(1:mp*16)); title('Unipolar NRZ');
subplot(4,2,2); pwelch(unrz,[],[],[],Fs,'power'); title('PSD of Unipolar NRZ');
% filter 0.046
unrz_f1 = conv(f1,unrz);
subplot(4,2,3); stem(unrz_f1(1:mp*16)); title('Unipolar NRZ fc = 0.046');
subplot(4,2,4); pwelch(unrz_f1,[],[],[],Fs,'power'); title('PSD of Unipolar NRZ fc = 0.046');
% filter 0.146
unrz_f2 = conv(f2,unrz);
subplot(4,2,5); stem(unrz_f2(1:mp*16)); title('Unipolar NRZ fc = 0.146');
subplot(4,2,6); pwelch(unrz_f2,[],[],[],Fs,'power'); title('PSD of Unipolar NRZ fc = 0.146');
% filter 0.246
unrz_f3 = conv(f3,unrz);
subplot(4,2,7); stem(unrz_f3(1:mp*16)); title('Unipolar NRZ fc = 0.246');
subplot(4,2,8); pwelch(unrz_f3,[],[],[],Fs,'power'); title('PSD of Unipolar NRZ fc = 0.246');
%% Filtre la señal Polar NRZ con los tres filtros. La salida del filtro representa la señal recibida
figure;
subplot(4,2,1); stem(pnrz(1:mp*16)); title('Polar NRZ');
subplot(4,2,2); pwelch(pnrz,[],[],[],Fs,'power'); title('PSD of Polar NRZ');
% filter 0.046
pnrz_f1 = conv(f1,pnrz);
subplot(4,2,3); stem(pnrz_f1(1:mp*16)); title('Polar NRZ fc = 0.046');
subplot(4,2,4); pwelch(pnrz_f1,[],[],[],Fs,'power'); title('PSD of Polar NRZ fc = 0.046');
% filter 0.146
pnrz_f2 = conv(f2,pnrz);
subplot(4,2,5); stem(pnrz_f2(1:mp*16)); title('Polar NRZ fc = 0.146');
subplot(4,2,6); pwelch(pnrz_f2,[],[],[],Fs,'power'); title('PSD of Polar NRZ fc = 0.046');
% filter 0.246
pnrz_f3 = conv(f3,pnrz);
subplot(4,2,7); stem(pnrz_f3(1:mp*16)); title('Polar NRZ fc = 0.246');
subplot(4,2,8); pwelch(pnrz_f3,[],[],[],Fs,'power'); title('PSD of Polar NRZ fc = 0.046');
%% Filtre la señal Manchester con los tres filtros. La salida del filtro representa la señal recibida
figure;
subplot(4,2,1); stem(xM(1:mp*16)); title('Manchester');
subplot(4,2,2); pwelch(xM,[],[],[],Fs,'power'); title('PSD of Manchester');
% filter 0.046
xM_f1 = conv(f1,xM);
subplot(4,2,3); stem(xM_f1(1:mp*16)); title('Manchester fc = 0.046');
subplot(4,2,4); pwelch(xM_f1,[],[],[],Fs,'power'); title('PSD of Manchester fc = 0.146');
% filter 0.146
xM_f2 = conv(f2,xM);
subplot(4,2,5); stem(xM_f2(1:mp*16)); title('Manchester fc = 0.146');
subplot(4,2,6); pwelch(xM_f2,[],[],[],Fs,'power'); title('PSD of Manchester fc = 0.146');
% filter 0.246
xM_f3 = conv(f3,xM);
subplot(4,2,7); stem(xM_f3(1:mp*16)); title('Manchester fc = 0.246');
subplot(4,2,8); pwelch(xM_f3,[],[],[],Fs,'power'); title('PSD of Manchester fc = 0.246');
% filter 0.8
xM_f5 = conv(f5,xM);
%% Filtre la señal Bipolar nRZ con los tres filtros. La salida del filtro representa la señal recibida
figure;
subplot(4,2,1); stem(bnrz(1:mp*16)); title('Bipolar NRZ');
subplot(4,2,2);pwelch(bnrz,[],[],[],Fs,'power'); title('PSD of Bipolar NRZ');
% filter 0.046
bnrz_f1 = conv(f1,bnrz);
subplot(4,2,3); stem(bnrz_f1(1:mp*16)); title('Bipolar NRZ fc = 0.046');
subplot(4,2,4);pwelch(bnrz_f1,[],[],[],Fs,'power'); title('PSD of Bipolar NRZ fc = 0.046');
% filter 0.146
bnrz_f2 = conv(f2,bnrz);
subplot(4,2,5); stem(bnrz_f2(1:mp*16)); title('Bipolar NRZ fc = 0.146');
subplot(4,2,6); pwelch(bnrz_f2,[],[],[],Fs,'power'); title('PSD of Bipolar NRZ fc = 0.146');
% filter 0.246
bnrz_f3 = conv(f3,bnrz);
subplot(4,2,7); stem(bnrz_f3(1:mp*16)); title('Bipolar NRZ fc = 0.246');
subplot(4,2,8); pwelch(bnrz_f3,[],[],[],Fs,'power'); title('PSD of Bipolar NRZ fc = 0.246');

%% instantes de observación & umbral de desicion

% % % % % % % % % % % % % % % % % % % % Unipolar NRZ 
y_s_unrz_f1 = unrz_f1(55:mp:end);
numel(y_s_unrz_f1) % 8718
y_s_unrz_f1 = y_s_unrz_f1(1:8712);
% scatterplot(y_s_unrz_f1)
sym_Rx_unrz_f1 = sign(y_s_unrz_f1 - 0.5);
bits_Rx_unrz_f1 = (sym_Rx_unrz_f1 +1)/2;

y_s_unrz_f2 = unrz_f2(55:mp:end);
numel(y_s_unrz_f2) % 8718
y_s_unrz_f2 = y_s_unrz_f2(1:8712);
% scatterplot(y_s_unrz_f2)
sym_Rx_unrz_f2 = sign(y_s_unrz_f2 - 0.5);
bits_Rx_unrz_f2 = (sym_Rx_unrz_f2 +1)/2;

y_s_unrz_f3 = unrz_f3(55:mp:end);
numel(y_s_unrz_f3) % 8718
y_s_unrz_f3 = y_s_unrz_f3(1:8712);
% scatterplot(y_s_unrz_f3)
sym_Rx_unrz_f3 = sign(y_s_unrz_f3 - 0.5);
bits_Rx_unrz_f3 = (sym_Rx_unrz_f3 +1)/2;

% % % % % % % % % % % % % % % % % % % % Polar NRZ 
y_s_pnrz_f1 = pnrz_f1(55:mp:end);
numel(y_s_pnrz_f1) % 8718
y_s_pnrz_f1 = y_s_pnrz_f1(1:8712);
% scatterplot(y_s_pnrz_f1)
sym_Rx_pnrz_f1 = sign(y_s_pnrz_f1);
bits_Rx_pnrz_f1 = (sym_Rx_pnrz_f1 +1)/2; % bits_Rx_pnrz_f1 es lo que se va a reconstruir

y_s_pnrz_f2 = pnrz_f2(55:mp:end);
numel(y_s_pnrz_f2) % 8718
y_s_pnrz_f2 = y_s_pnrz_f2(1:8712);
% scatterplot(y_s_pnrz_f2)
sym_Rx_pnrz_f2 = sign(y_s_pnrz_f2);
bits_Rx_pnrz_f2 = (sym_Rx_pnrz_f2 +1)/2; % bits_Rx_pnrz_f2 es lo que se va a reconstruir

y_s_pnrz_f3 = pnrz_f3(55:mp:end);
numel(y_s_pnrz_f3) % 8718
y_s_pnrz_f3 = y_s_pnrz_f3(1:8712);
% scatterplot(y_s_pnrz_f3)
sym_Rx_pnrz_f3 = sign(y_s_pnrz_f3);
bits_Rx_pnrz_f3 = (sym_Rx_pnrz_f3 +1)/2; % bits_Rx_pnrz_f3 es lo que se va a reconstruir

% % % % % % % % % % % % % % % % % % % % Manchester 
y_s_xM_f1 = xM_f1(53:mp:end);
numel(y_s_xM_f1) % 8718
y_s_xM_f1 = y_s_xM_f1(1:8712);
% scatterplot(y_s_xM_f1)
sym_Rx_xM_f1 = sign(y_s_xM_f1);
bits_Rx_xM_f1 = (sym_Rx_xM_f1 +1)/2;
% not en el receptor de los bits debido a que escogi una convención de Manchester opcion B 
bits_Rx_xM_f1 = not(bits_Rx_xM_f1)

y_s_xM_f2 = xM_f2(53:mp:end);
numel(y_s_xM_f2) % 8718
y_s_xM_f2 = y_s_xM_f2(1:8712);
% scatterplot(y_s_xM_f2)
sym_Rx_xM_f2 = sign(y_s_xM_f2);
bits_Rx_xM_f2 = (sym_Rx_xM_f2 +1)/2; 
% not en el receptor de los bits debido a que escogi una convención de Manchester opcion B 
bits_Rx_xM_f2 = not(bits_Rx_xM_f2)

y_s_xM_f3 = xM_f3(53:mp:end);
numel(y_s_xM_f3) % 8718
y_s_xM_f3 = y_s_xM_f3(1:8712);
% scatterplot(y_s_xM_f3)
sym_Rx_xM_f3 = sign(y_s_xM_f3);
bits_Rx_xM_f3 = (sym_Rx_xM_f3 +1)/2;
% not en el receptor de los bits debido a que escogi una convención de Manchester opcion B 
bits_Rx_xM_f3 = not(bits_Rx_xM_f3) 

%% 
y_s_xM_f5 = xM_f5(53:mp:end);
numel(y_s_xM_f5) % 8718
y_s_xM_f5 = y_s_xM_f5(1:8712);
% scatterplot(y_s_xM_f5)
sym_Rx_xM_f5 = sign(y_s_xM_f5);
bits_Rx_xM_f5 = (sym_Rx_xM_f5 +1)/2;
% not en el receptor de los bits debido a que escogi una convención de Manchester opcion B 
bits_Rx_xM_f5 = not(bits_Rx_xM_f5)

%%
% % % % % % % % % % % % % % % % % % % % Bipolar NRZ 
y_s_bnrz_f1 = bnrz_f1(55:mp:end);
numel(y_s_bnrz_f1) % 8718
y_s_bnrz_f1 = y_s_bnrz_f1(1:8712);
% scatterplot(y_s_bnrz_f1)
y_s_bnrz_f1 = abs(y_s_bnrz_f1); 
sym_Rx_bnrz_f1 = sign(y_s_bnrz_f1 - 0.5);
bits_Rx_bnrz_f1 = (sym_Rx_bnrz_f1 +1)/2;

y_s_bnrz_f2 = bnrz_f2(55:mp:end);
numel(y_s_bnrz_f2) % 8718
y_s_bnrz_f2 = y_s_bnrz_f2(1:8712);
% scatterplot(y_s_bnrz_f2)
y_s_bnrz_f2 = abs(y_s_bnrz_f2);  
sym_Rx_bnrz_f2 = sign(y_s_bnrz_f2 - 0.5);
bits_Rx_bnrz_f2 = (sym_Rx_bnrz_f2 +1)/2;

y_s_bnrz_f3 = bnrz_f3(55:mp:end);
numel(y_s_bnrz_f3) % 8718
y_s_bnrz_f3 = y_s_bnrz_f3(1:8712);
% scatterplot(y_s_bnrz_f3)
y_s_bnrz_f3 = abs(y_s_bnrz_f3); 
sym_Rx_bnrz_f3 = sign(y_s_bnrz_f3 - 0.5);
bits_Rx_bnrz_f3 = (sym_Rx_bnrz_f3 +1)/2;
%%
% bR = reshape(bits_Rx_unrz_f1,[8,1089]);
% % isequal(bR,b)
% bR = bR';
% lena = reshape(bi2de(bR),size(lenarec));
% isequal(lenarec,lena)
% figure;
% imshow(uint8(lena))
% %% 
% %%
% bR = reshape(bits_Rx_unrz_f2,[8,1089]);
% % isequal(bR,b)
% bR = bR';
% lena = reshape(bi2de(bR),size(lenarec));
% isequal(lenarec,lena)
% figure;
% imshow(uint8(lena))
%% Para reconstruir la imagen hay que cambiar el primer argumento del reshape por la variable
%%
bR = reshape(bits_Rx_xM_f5,[8,1089]);
% isequal(bR,b)
bR = bR';
lena = reshape(bi2de(bR),size(lenarec));
isequal(lenarec,lena)
figure;
imshow(uint8(lena))

% bits_Rx_bnrz_f3
%% Osciloscopio
% Fs = 44100;
% Fs = 48000;
Fs = 96000;
soundsc(xM, Fs)
% The bit rate is Rb = Rs =  Fs / mp
Rb = Rs
%% ejercicio 3
pulso = triang(100); % pulso triangular de 100 muestras
Ep = sum(pulso.*pulso) % energía del pulso
E1 = max(conv(pulso,fliplr(pulso))) % máximo de la convolución
E2 = conv(pulso,fliplr(pulso));
E2(100) % elemento 100 de la convolución