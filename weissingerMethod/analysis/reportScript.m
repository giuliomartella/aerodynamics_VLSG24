clear; clc; close all
% Generate report plots

addpath('../'); % load functions
[wingC, tailC] = cessna(); % load data
[wingG, tailG] = gliderK8();



%% Build Elements
plotFlag = false;
wingC = buildElement(wingC, plotFlag);
tailC = buildElement(tailC, plotFlag);
wingG = buildElement(wingG, plotFlag);
tailG = buildElement(tailG, plotFlag);

%% Compute
alpha_range = deg2rad(-5:1:12);
cL = zeros(length(alpha_range), 2);
cD = cL;
for i = 1:length(alpha_range)
    alpha = alpha_range(i);
    uInf = [cos(alpha), -sin(alpha), 0; sin(alpha), cos(alpha), 0; 0, 0, 1] * [1.0; 0.0; 0.0];

    wingC = buildLinearSystem(uInf, wingC);
    wingC = postprocessing(wingC);
    wingG = buildLinearSystem(uInf, wingG);
    wingG = postprocessing(wingG);

    cL(i, 1) = wingC.cL;
    cD(i, 1) = wingC.cD;
    cL(i, 2) = wingG.cL;
    cD(i, 2) = wingG.cD;

    % Add Tail

    [wingC, tailC] = buildLinearSystem(uInf, wingC, tailC);
    wingC = postprocessing(wingC, tailC);
    [wingG, tailG] = buildLinearSystem(uInf, wingG, tailC);
    wingG = postprocessing(wingG, tailG);

    cL(i, 3) = wingC.cL;
    cD(i, 3) = wingC.cD;
    cL(i, 4) = wingG.cL;
    cD(i, 4) = wingG.cD;
end

%% Polar
cdElliptic = @(cl, AR) cl.^2 / pi / AR; 
figure(6);
plot(cD(:, 1), cL(:, 1), 'bo-', 'LineWidth', 1.5, 'MarkerSize', 6);
grid on
hold on
plot(cD(:, 2), cL(:, 2), 'ro-', 'LineWidth', 1.5, 'MarkerSize', 6);
plot(cdElliptic(linspace(-0.5, 1.3), wingC.AR), linspace(-0.5, 1.3), 'b', 'LineWidth', 0.5);
plot(cdElliptic(linspace(-0.5, 1.3), wingG.AR), linspace(-0.5, 1.3),'r', 'LineWidth', 0.5);
plot(cD(:,3), cL(:,3), '--b^', 'LineWidth', 1.5, 'MarkerSize', 6, 'MarkerFaceColor', 'blue'); % Cessna con coda
plot(cD(:,4), cL(:,4), '--r^', 'LineWidth', 1.5, 'MarkerSize', 6, 'MarkerFaceColor', 'red'); % Glider con coda


% Set labels and title with LaTeX formatting
xlabel('$C_D$ (Induced Drag Coefficient)', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('$C_L$ (Lift Coefficient)', 'Interpreter', 'latex', 'FontSize', 14);
title('Lift-Drag Polar Curve', 'Interpreter', 'latex', 'FontSize', 16);

% Adjust axis for better visualization
axis tight;

% Optional: Add a legend if needed
legend('Cessna wing', 'Glider wing', 'Elliptic distribution with Cessna AR ', 'Elliptic distribution with Glider AR ','Cessna with tail', 'Glider with tail', 'Interpreter', 'latex', 'FontSize', 12, 'Location', 'best');

% saveas(gcf, 'polar.pdf');

 %% Cl Alpha
figure(7);
plot(rad2deg(alpha_range), cL(:,1), '-bo', 'LineWidth', 1.5, 'MarkerSize', 6);
hold on
plot(rad2deg(alpha_range), cL(:,2), '-ro', 'LineWidth', 1.5, 'MarkerSize', 6);
plot(rad2deg(alpha_range), alpha_range*2*pi, 'g', 'LineWidth', 0.5);
grid on
plot(rad2deg(alpha_range), cL(:,3), '--b^', 'LineWidth', 1.5, 'MarkerSize', 6, 'MarkerFaceColor', 'blue'); % Cessna con coda
plot(rad2deg(alpha_range), cL(:,4), '--r^', 'LineWidth', 1.5, 'MarkerSize', 6, 'MarkerFaceColor', 'red'); % Glider con coda


% Set labels and title with LaTeX formatting
xlabel('$Alpha$ (AoA)', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('$C_L$ (Lift Coefficient)', 'Interpreter', 'latex', 'FontSize', 14);
title('Lift-AoA Polar Curve', 'Interpreter', 'latex', 'FontSize', 16);

% Adjust axis for better visualization
axis tight;

% Add a legend
legend('Cessna wing','Glider wing', 'Two Dimensional Thin Airfoil','Cessna with tail', 'Glider with tail', 'Interpreter', 'latex', 'FontSize', 12, 'Location', 'best');

% saveas(gcf, 'cl_alpha.pdf');

%% Circulation distribution, AoA = 3°

alpha = deg2rad(3.0);
uInf = [cos(alpha), -sin(alpha), 0; sin(alpha), cos(alpha), 0; 0, 0, 1] * [1.0; 0.0; 0.0];
wingC = buildLinearSystem(uInf, wingC, tailC);
wingC = postprocessing(wingC);
wingG = buildLinearSystem(uInf, wingG, tailG);
wingG = postprocessing(wingG);

N = length(wingC.gammaDistribution); % Sine functions number
S = zeros(length(wingC.s), N);
for n = 1:N
    S(:, n) = sin(n * pi * (wingC.s + 0.5));
end
a = S \ wingC.gammaDistribution'; % Fit Coefficients
% Build series
gammaEllipticC = S(:,1) * a(1);

N = length(wingG.gammaDistribution); % Sine functions number
S = zeros(length(wingG.s), N);
for n = 1:N
    S(:, n) = sin(n * pi * (wingG.s + 0.5));
end
a = S \ wingG.gammaDistribution'; % Fit Coefficients
% Build series
gammaEllipticG = S(:,1) * a(1);

figure(8);
plot(wingC.s, wingC.gammaDistribution, '-ro', 'LineWidth', 1.5, 'MarkerSize', 6);
grid on;
hold on
plot(wingG.s, wingG.gammaDistribution, '-bo', 'LineWidth', 1.5, 'MarkerSize', 6);
plot(wingC.s, gammaEllipticC, '-r', 'LineWidth', 0.5, 'MarkerSize', 6);
plot(wingG.s, gammaEllipticG, '-b', 'LineWidth', 0.5, 'MarkerSize', 6);
scatter([-0.5; 0.5], [0; 0]);


% Set labels and title with LaTeX formatting
xlabel('$s$ (Span Fraction)', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('$\Gamma$ (Circulation)', 'Interpreter', 'latex', 'FontSize', 14);
title('Circulation distribution', 'Interpreter', 'latex', 'FontSize', 16);

% Adjust axis for better visualization
axis tight;

% Optional: Add a legend if needed
legend('Cessna','Glider', 'Cessna Elliptic distribution','Glider Elliptic distribution',  'Interpreter', 'latex', 'FontSize', 12, 'Location', 'best');

% saveas(gcf, 'circulation.pdf');

