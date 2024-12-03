close all
clear
clc

% Loading struct with all the necessary data
% Points : x and y coordinates of  the profile
% Angles : angles for which cp has been retrieved
% Cp : matrix with cp data for all angles. Each column is a different AOA


load("cpData.mat");

% Finding leading edge
ba = find(~cpData.Points(:,1));
ba = ba(1);

% Cycle analyzing cp minimum and leading edge difference
cpPeak = zeros(length(cpData.Angles),1);
pos = zeros(length(cpData.Angles),1);
leCpDiff = zeros(length(cpData.Angles),1);

for ii = 1:length(cpData.Angles)
    cpup = interp1(cpData.Points(1:ba,1), cpData.Cp(1:ba,ii),linspace(0,1,1000));
    cpdown = interp1(cpData.Points(ba+1:end,1), cpData.Cp(ba+1:end,ii),linspace(0,1,1000));
    cpDiff = cpdown - cpup;
    [cpPeak(ii),pos(ii)] = max(cpDiff);
    leCpDiff(ii) = abs( cpData.Cp(ba-1,ii) - cpData.Cp(ba+1,ii));
end


% Finding minimum peak
[minCpPeak,I] = min(cpPeak);
alphaTheodorsen = cpData.Angles(I);
fprintf("Theodorsen angle according to minimum cp difference: %.4f°", cpData.Angles(I))
fprintf("\n Minimum peak is of %.4f and is achieved at %.4f m along the chord\n\n", minCpPeak, pos(I)/1000)

[~,ind] = min(leCpDiff);
fprintf("Theodorsen angle according to minimum cp difference on leading edge: %.4f°", cpData.Angles(ind))

%% CP PLOT - Choose an angle
Angle = alphaTheodorsen;


jj = find (~(cpData.Angles - Angle));

figure;
hold on
axis equal
grid on
plot(cpData.Points(1:ba,1),cpData.Points(1:ba,2),'b')
plot(cpData.Points(ba+1:end,1),cpData.Points(ba+1:end,2),'r')

cpup = interp1(cpData.Points(1:ba,1), cpData.Cp(1:ba,jj),linspace(0,1,1000));
cpdown = interp1(cpData.Points(ba+1:end,1), cpData.Cp(ba+1:end,jj),linspace(0,1,1000));
cpDiff = cpdown - cpup;

plot(linspace(0,1,1000), cpDiff)

plot(cpData.Points(1:ba,1),cpData.Cp(1:ba,jj),'b')
plot(cpData.Points(ba+1:end,1),cpData.Cp(ba+1:end,jj),'r')

legend("Dorso","Ventre", "cpDiff", "cpDorso","cpVentre")
