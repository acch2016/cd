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
bits=b(:);

mp = 10; % samples per pulse
Fs = 96000;
Ts = 1/Fs;
% The bit rate is Rb= Rs= Fs / mp, because 1 bit= 1 %symbol and % every symbol has mpsamples per bit
Rs = Fs / mp;

p = input('elegir pulso base (1 cuadrado / 2 medio-seno RZ): ');
switch p
    case 1
        pbase = rectwin(mp); % Pulso completo UNRZ PNRZ BNRZ;
        pm = [ones(1,mp/2) -ones(1,mp/2)]; %pulso base Manchester opcion A
    case 2
        mb2 = round(mp/2);
        n = 0:mb2-1; w0=pi/(mb2);
        half_sin = sin(w0*n); % Probar con medio-seno RZ como en Zigbee
        pm = [half_sin -half_sin]; pbase = pm;
end

linecode = input('elegir lineCode (unrz/pnrz/bnrz/m): ','s');
if linecode == 'unrz'
    
    wvtool(pbase);
    s = zeros(1,numel(bits)*mp);
    s(1:mp:end) = bits;
    subplot(3,1,1); stem(s(1:mp*16)); title('tren de impulsos');
    unrz = conv(pbase,s);
    subplot(3,1,2); plot(unrz(1:mp*16)); title('Unipolar NRZ');
    subplot(3,1,3); pwelch(unrz,[],[],[],Fs,'power'); title('PSD of Unipolar NRZ');
    
    % filtrar
    figure; subplot(2,2,1); stem(unrz(1:mp*16)); title('Unipolar NRZ');
    subplot(2,2,2); pwelch(unrz,[],[],[],Fs,'power'); title('PSD of Unipolar NRZ');
    rx = conv(fil,unrz); % señal recibida
    subplot(2,2,3); stem(rx(1:mp*16)); title(['Unipolar NRZ fc = ',num2str(fc)]);
    subplot(2,2,4); pwelch(rx,[],[],[],Fs,'power'); title(['PSD of Unipolar NRZ fc = ',num2str(fc)]);
    
    % % % % % % % % % % % % % % % instantes de observación & umbral de desicion
    y_s = rx(55:mp:end);
%     numel(y_s) % 8718
    y_s = y_s(1:8712);
    % scatterplot(y_s_unrz_f1)
    sym_Rx = sign(y_s - 0.5);
    bits_Rx = (sym_Rx +1)/2;

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
if linecode == 'bnrz'
    
    wvtool(pbase);
    am = mod(1:length(bits(bits == 1)), 2); % 0s y 1s del tamaño de los bits en 1 en la variable bits
    am(am == 0) = -1; % alternate mark (vector de -1s y 1s)
    bits1 = bits; % para no modificar la variable original que contine los bits
    bits1(bits == 1) = am; % alternate mark inversion AMI (tres estados)
    s = zeros(1, numel(bits1)*mp);
    s(1:mp:end) = bits1;
    subplot(3,1,1); stem(s(1:mp*16)); title('tren de impulsos');
    bnrz = conv(pbase,s);
    subplot(3,1,2); plot(bnrz(1:mp*16)); title('Bipolar NRZ');
    subplot(3,1,3); pwelch(bnrz,[],[],[],Fs,'power'); title('PSD of Bipolar NRZ');
    % figure;powerbw(bnrz,Fs) % para obtener el bw de -3dB
    
    figure;
    subplot(2,2,1); stem(bnrz(1:mp*16)); title('Bipolar NRZ');
    subplot(2,2,2); pwelch(bnrz,[],[],[],Fs,'power'); title('PSD of Bipolar NRZ');
    rx = conv(fil,bnrz); % señal recibida
    subplot(2,2,3); stem(rx(1:mp*16)); title(['Bipolar NRZ fc = ', num2str(fc)]);
    subplot(2,2,4); pwelch(rx,[],[],[],Fs,'power'); title(['PSD of Bipolar NRZ fc = ', num2str(fc)]);
    
    % % % % % % % % % % % % % % % instantes de observación & umbral de desicion
    y_s = rx(55:mp:end);
    numel(y_s) % 8718
    y_s = y_s(1:8712);
    % scatterplot(y_s_bnrz_f1)
    y_s = abs(y_s);
    sym_Rx = sign(y_s - 0.5);
    bits_Rx = (sym_Rx +1)/2;

end
if linecode == 'm'
    
    wvtool(pm);
    sym = (bits * 2)-1;
    s = zeros(1,numel(bits)*mp);% tren de impulsos
    s(1:mp:end) = sym; %tren
    subplot(3,1,1); stem(s(1:mp*16)); title('tren de impulsos');
    xM = conv(pm,s);
    %     subplot(3,1,3); pwelch(xM,200,100,200,Fs,'power'); title('PSD of Manchester');
    subplot(3,1,3); pwelch(xM,[],[],[],Fs,'power'); title('PSD of Manchester');
    subplot(3,1,2); plot(xM(1:mp*16));title('Manchester');
    
    figure;
    subplot(2,2,1); stem(xM(1:mp*16)); title('Manchester');
    subplot(2,2,2); pwelch(xM,[],[],[],Fs,'power'); title('PSD of Manchester');
    rx = conv(fil,xM); % señal recibida
    subplot(2,2,3); stem(rx(1:mp*16)); title(['Manchester fc = ', num2str(fc)]);
    subplot(2,2,4); pwelch(rx,[],[],[],Fs,'power'); title(['PSD of Manchester fc = ', num2str(fc)]);
    
    % % % % % % % % % % % % % % % instantes de observación & umbral de desicion
    y_s = rx(53:mp:end);
%     numel(y_s) % 8718
    y_s = y_s(1:8712);
    % scatterplot(y_s_xM_f1)
    sym_Rx = sign(y_s);
    bits_Rx = (sym_Rx +1)/2;
    % not en el receptor de los bits no se necestita debido a que escogio convención Manchester opcion A
    %     bits_Rx = not(bits_Rx);

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


