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



dyf = @(x) spline(meanLine.xMl, meanLine.dy, x);

f_th = @(t) dyf(t)./sqrt(0.25 - t.^2);  
aTh = 1/pi *integral(f_th,meanLine.xMl(1), meanLine.xMl(end));  % using integral results are kinda better
aTh = rad2deg(aTh)

f_a0 = @(t) dyf(t).* sqrt((0.5 + t) / (0.5 -t));  
a0 = 2/pi *integral(f_a0,meanLine.xMl(1), meanLine.xMl(end));
a0 = rad2deg(a0)


