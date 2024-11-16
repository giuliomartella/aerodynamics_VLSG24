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
plot(meanLine(:,1), meanLine(:,2))


h = 1e-10; % PROBLEM!! the integral result depens only on this element
dy = meanLine(:,3);
t = linspace(-0.5 +h, 0.5 -h, numel(dy))';
aTh = trapz(t, dy./sqrt(0.25 - t.^2)) / pi
aTh = rad2deg(aTh)



dy = meanLine(:,3);
t = linspace(-0.5, 0.5, numel(dy))';
dyf = @(x) interp1(t, dy, x);
f = @(t) dyf(t)./sqrt(0.25 - t.^2);  

aTh = 1/pi *integral(f,-0.5 ,0.5)  % using integral results are kinda better
aTh = rad2deg(aTh)