close all
clear
clc

% Loading struct with all the necessary data
% Points : x and y coordinates of  the profile
% Angles : angles for which cp has been retrieved
% Cp : matrix with cp data for all angles. Each column is a different AOA


%load("alfa_cp/cpData.mat");
if ismac
    load("alfa_cp\cpData.mat");
else
    load("alfa_cp/cpData.mat");
end

% Finding indexes of leading edge and trailing edge
ba = find(~cpData.Points(:,1));
if length(ba) > 1
    baUp = ba(1) ;
    baDown = ba(2) ;
else 
    baUp = ba(1) - 1;
    baDown = ba(1) + 1;
end
% Cycle analyzing cp minimum and leading edge difference
cpPeak = zeros(length(cpData.Angles),1);
pos = zeros(length(cpData.Angles),1);
leCpDiff = zeros(length(cpData.Angles),1);

for ii = 1:length(cpData.Angles)
    cpUp = interp1(cpData.Points(1:baUp,1), cpData.Cp(1:baUp,ii),linspace(0,1,1000));
    cpDown = interp1(cpData.Points(baDown:end,1), cpData.Cp(baDown:end,ii),linspace(0,1,1000));
    cpDiff = cpDown - cpUp;
    [cpPeak(ii),pos(ii)] = max(abs(cpDiff));
    leCpDiff(ii) = abs( cpData.Cp(baDown,ii) - cpData.Cp(baUp,ii));
end


% Finding minimum peak
[minCpPeak,IndexTheo1] = min(cpPeak);

fprintf("Theodorsen angle according to minimum cp difference: %.4f°", cpData.Angles(IndexTheo1))
fprintf("\n Minimum peak is of %.4f and is achieved at %.4f m along the chord\n\n", minCpPeak, pos(IndexTheo1)/1000)

[~,IndexTheo2] = min(leCpDiff);
fprintf("Theodorsen angle according to minimum cp difference on leading edge: %.4f° \n\n", cpData.Angles(IndexTheo2))

%% CP PLOT - Minimized cp peak 
Angle = cpData.Angles(IndexTheo1);

jj = find (~(cpData.Angles - Angle));

figure;
hold on
axis equal
grid on
plot(cpData.Points(1:baUp,1),cpData.Points(1:baUp,2),'b')
plot(cpData.Points(baDown:end,1),cpData.Points(baDown:end,2),'r')

cpUp = interp1(cpData.Points(1:baUp,1), cpData.Cp(1:baUp,jj),linspace(0,1,1000));
cpDown = interp1(cpData.Points(baDown:end,1), cpData.Cp(baDown:end,jj),linspace(0,1,1000));
cpDiff = cpDown - cpUp;

plot(linspace(0,1,1000), cpDiff)

plot(cpData.Points(1:baUp,1),cpData.Cp(1:baUp,jj),'b')
plot(cpData.Points(baDown:end,1),cpData.Cp(baDown:end,jj),'r')

legend("Dorso","Ventre", "cpDiff", "cpDorso","cpVentre")
title("Minimized cp Peak")

%% CP PLOT - Minimized leading edge cp difference
Angle = cpData.Angles(IndexTheo2);

jj = find (~(cpData.Angles - Angle));

figure;
hold on
axis equal
grid on
plot(cpData.Points(1:baUp,1),cpData.Points(1:baUp,2),'b')
plot(cpData.Points(baDown:end,1),cpData.Points(baDown:end,2),'r')

cpUp = interp1(cpData.Points(1:baUp,1), cpData.Cp(1:baUp,jj),linspace(0,1,1000));
cpDown = interp1(cpData.Points(baDown:end,1), cpData.Cp(baDown:end,jj),linspace(0,1,1000));
cpDiff = cpDown - cpUp;

plot(linspace(0,1,1000), cpDiff)

plot(cpData.Points(1:baUp,1),cpData.Cp(1:baUp,jj),'b')
plot(cpData.Points(baDown:end,1),cpData.Cp(baDown:end,jj),'r')

legend("Dorso","Ventre", "cpDiff", "cpDorso","cpVentre")
title("Minimized cp on le")