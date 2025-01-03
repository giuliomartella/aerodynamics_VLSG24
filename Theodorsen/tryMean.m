clear
clc
close all

% Add the directory containing airfoil data to the path
addpath("mat_functions/");

% Import airfoil profile data
filename = 'k2.txt';
airfoil = importXfoilProfile(filename);

% Compute mean camber line and its derivative
meanLine = meanLineFun(airfoil);

% Plot the airfoil, mean camber line, and its derivative
figure
plot(meanLine.x, meanLine.y) % Original airfoil
hold on
axis equal
grid on
plot(meanLine.xMl, meanLine.yMl, 'r') % Mean camber line
plot(meanLine.xMl, meanLine.dy, 'g') % Camber line derivative

% Spline interpolation for camber line derivative
dyf = @(x) spline(meanLine.xMl, meanLine.dy, x);

%  Note: Rotate the profile back by the angle tau to return it to the 
%  original reference frame used in the Xfoil .dat file.

% Compute Theodorsen's angle of attack
f_th = @(t) dyf(t) ./ sqrt(0.25 - t.^2);  
aTh = 1/pi * integral(f_th, -0.5, 0.5) + meanLine.tau;
aThDeg = rad2deg(aTh);
disp(['Theodorsen angle of attack (degrees): ', num2str(aThDeg)])

% Compute zero-lift angle of attack
f_a0 = @(t) sqrt(0.5 + t) .* dyf(t) ./ sqrt(0.5 - t);
a0 = 2/pi * integral(f_a0, -0.5, 0.5) + meanLine.tau;
a0Deg = rad2deg(a0);
disp(['Zero-lift angle of attack (degrees): ', num2str(a0Deg)])



%  Integrals computed through theta variable

% f_th = @(t) dyf(-0.5 * cos(t));  
% aTh = 1/pi *integral(f_th, 0, pi);  
% aThDeg = rad2deg(aTh)
% 
% f_a0 = @(t) dyf(-0.5 * cos(t)).* cos(t);  
% a0 = aTh - 1/pi *integral(f_a0, 0, pi) + meanLine.tau;
% a0Deg = rad2deg(a0)
% 
% f_a0 = @(t) dyf(-0.5 * cos(t)).* (1-cos(t));  
% a0 = 1/pi *integral(f_a0, 0, pi) + meanLine.tau;
% a0Deg = rad2deg(a0)


