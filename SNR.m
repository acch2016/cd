clc;
clear all;
x = rand(1, 1e6);%renglones, columnas
min(x)
hist(x, 10)
%promedio es 0.5
%plot()
promedio_manual = sum(x)/numel(x)
promedio = mean(x)
x1 = x-promedio;
mean(x1)
min(x1)
figure;hist(x1, 10)
var(x)
sum((x-mean(x)).^2)/numel(x) % misma operacion para todos los elementos
xg = randn(1, 1e6);%n de normal
hist(xg, 100)
min(xg)
max(xg)
mean(xg)
var(xg)
Pxg =sum((xg-mean(xg)).^2)/numel(xg)