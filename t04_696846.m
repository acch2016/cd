% mp = 9; % muestras por pulso
% pbase = triang(mp); % pulso base triangular de mp muestras
% t = 0:8;
% b = [1 0 1 0 1 1 0 0];
% stem(t,pbase)
% ylim([0 1])
% s = [1 0 0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 1];
% s1 = b;
% s1(s1==0) = -1;
% s = zeros(1,(numel(b)-1)*mp+1);
% s(1:mp:end) = s1;
% x = conv(s,pbase);
% stem(x)
load lena512.mat
whos('-file','lena512.mat')
imshow(uint8(lena512))
lenarec=lena512(252:284,318:350);
imshow(uint8(lenarec))
b=de2bi(lenarec,8);
b=b';
bits=b(:);