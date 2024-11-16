clear
clc
close all

addpath airfoils
filename = 'k2.dat';
airfoil = importXfoilProfile(filename);

meanLine = meanLineFun(airfoil);

figure
scatter(airfoil.x, airfoil.y)
hold on
axis equal
grid on
plot(meanLine(:,1), meanLine(:,2))