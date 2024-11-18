clear
clc
close all

addpath airfoils
filename = 'k2_1.dat';
airfoil = importXfoilProfile(filename);

meanLine = meanLineFun(airfoil);

figure
plot(airfoil.x, airfoil.y)
hold on
axis equal
grid on
plot(meanLine.xMl, meanLine.yMl)
plot(meanLine.xMl, meanLine.dy)



t = linspace(-0.5, 0.5, numel(meanLine.dy))';
dyf = @(x) spline(t, meanLine.dy, x);

f_th = @(t) dyf(t)./sqrt(0.25 - t.^2);  

aTh = 1/pi *integral(f_th,-0.5 ,0.5);  % using integral results are kinda better
aTh = rad2deg(aTh)

f_a0 = @(t) dyf(t).* sqrt((0.5 + t) / (0.5 -t));  
a0 = 2/pi *integral(f_a0,-0.5 ,0.5);
a0 = rad2deg(a0)


