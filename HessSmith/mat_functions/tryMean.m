clear
clc
close all

addpath airfoils
filename = 'k2.dat';
airfoil = importXfoilProfile(filename);

meanLine = meanLineFun(airfoil);

figure
plot(meanLine.x, meanLine.y)
hold on
axis equal
grid on
plot(meanLine.xMl, meanLine.yMl)
plot(meanLine.xMl, meanLine.dy)


h = 1e-6;
t = linspace(-0.5+h, 0.5-h, numel(meanLine.dy))';
dyf = @(x) interp1(t, meanLine.dy, x);

f_th = @(t) dyf(t)./sqrt(0.25 - t.^2);  
aTh = 1/pi *integral(f_th,t(1), t(end));  % using integral results are kinda better
aTh = rad2deg(aTh)

f_a0 = @(t) dyf(t).* sqrt((0.5 + t) / (0.5 -t));  
a0 = 2/pi *integral(f_a0,t(1), t(end));
a0 = rad2deg(a0)


