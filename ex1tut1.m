clear all;close all;
% señal muy similar a la analógica
Ts = 1/1000;
t = 0:Ts:1;
x = cos(2*pi*t*10);
% señal con muestreo natural 19 veces más lento (53Hz)
tm = t(1:19:end);
xm = cos(2*pi*tm*10);
% señal muestreada con sample-and-hold
ts = 19;
xs = zeros(1,length(t));
for i=1:length(t)
if( rem(i,ts)==1 )
tmp = x(i);
end
xs(i) = tmp;
end
% cuantificación
M = 64;
int = (max(xs)-min(xs))/M;
m = (min(xs)+int/2):int:(max(xs)-int/2);
xq = zeros(1,length(t));
for i=1:length(t)
[tmp k] = min(abs(xs(i)-m));
xq(i) = m(k);
end
% diferencia
xd = xs - xq;
% graficas
figure(1)
plot(t,x)  %analog signal
hold on
stem(tm,xm,'r','filled') %natural sampling
figure(2)
plot(t,xs)  %sample and hold
figure(3)
plot(t,xq) %quantified
figure(4)
plot(t,xd) %difference

