clear all; clc;
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
% b=de2bi(lenarec,8,'left-msb'); % % % % % % % version profesor
b=b';
bits=b(:);
linecode = input('elegir lineCode (unrz/pnrz/bnrz/m): ','s');

mp = 10; % samples per pulse
Fs = 96000;
Ts = 1/Fs;
% The bit rate is Rb= Rs= Fs / mp, because 1 bit= 1 %symbol and % every symbol has mpsamples per bit
Rs = Fs / mp;
if linecode == 'm'
%     mp = 10;
    mb2 = round(mp/2);
    n = 0:mb2-1; w0=pi/(mb2);
    half_sin = sin(w0*n); % Probar con medio-seno RZ como en Zigbee
    pm = [half_sin -half_sin];
    % pM = [ones(1,mp/2) -ones(1,mp/2)]; %pulso base
    wvtool(pm);
    sym = (bits * 2)-1;
    s = zeros(1,numel(bits)*mp);% tren de impulsos
    s(1:mp:end) = sym; %tren
    subplot(3,1,1); stem(s(1:mp*16)); title('tren de impulsos');
    xM = conv(pm,s);
    subplot(3,1,3); pwelch(xM,200,100,200,Fs,'power'); title('PSD of Manchester');
    subplot(3,1,3); pwelch(xM,[],[],[],Fs,'power'); title('PSD of Manchester');
    subplot(3,1,2); plot(xM(1:mp*16));title('Manchester');
end




