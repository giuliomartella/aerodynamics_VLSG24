clear; clc; close all
% Generate report plots

addpath('../'); % load functions
[wingC, tailC] = cessna(); % load data
[wingG, tailG] = gliderK8();



% Build Elements
plotFlag = false;
wingC = buildElement(wingC, true);
tailC = buildElement(tailC, plotFlag);
wingG = buildElement(wingG, false);
tailG = buildElement(tailG, plotFlag);

%% Compute
alpha_range = deg2rad(-4:2:15);
cL = zeros(length(alpha_range), 6);
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

    % Add Tail 0° angle

    [wingC, tailC] = buildLinearSystem(uInf, wingC, tailC);
    wingC = postprocessing(wingC, tailC);
    [wingG, tailG] = buildLinearSystem(uInf, wingG, tailG);
    wingG = postprocessing(wingG, tailG);

    cL(i, 3) = wingC.cL;
    cD(i, 3) = wingC.cD;
    cL(i, 4) = wingG.cL;
    cD(i, 4) = wingG.cD;
end

% Add Tail -2° angle
tailC.twistZero = deg2rad(-2.0);
tailG.twistZero = deg2rad(-2.0);
tailC = buildElement(tailC, true);
tailG = buildElement(tailG, false);

for i = 1:length(alpha_range)
    alpha = alpha_range(i);
    uInf = [cos(alpha), -sin(alpha), 0; sin(alpha), cos(alpha), 0; 0, 0, 1] * [1.0; 0.0; 0.0];

    [wingC, tailC] = buildLinearSystem(uInf, wingC, tailC);
    wingC = postprocessing(wingC, tailC);
    [wingG, tailG] = buildLinearSystem(uInf, wingG, tailG);
    wingG = postprocessing(wingG, tailG);

    cL(i, 5) = wingC.cL;
    cD(i, 5) = wingC.cD;
    cL(i, 6) = wingG.cL;
    cD(i, 6) = wingG.cD;
end


%% Polar
cdElliptic = @(cl, AR) cl.^2 / pi / AR; 
figure();
plot(cD(:, 1), cL(:, 1), 'bo-', 'LineWidth', 1.5, 'MarkerSize', 6); % cessna
grid on
hold on
plot(cD(:, 2), cL(:, 2), 'ro-', 'LineWidth', 1.5, 'MarkerSize', 6); % glider
plot(cD(:,3), cL(:,3), '--b^', 'LineWidth', 1.5, 'MarkerSize', 6, 'MarkerFaceColor', 'blue'); % Cessna 0
plot(cD(:,4), cL(:,4), '--r^', 'LineWidth', 1.5, 'MarkerSize', 6, 'MarkerFaceColor', 'red'); % Glider 0
plot(cD(:,5), cL(:,5), '--g^', 'LineWidth', 1.5, 'MarkerSize', 6, 'MarkerFaceColor', 'blue'); % Cessna -5
plot(cD(:,6), cL(:,6), '--m^', 'LineWidth', 1.5, 'MarkerSize', 6, 'MarkerFaceColor', 'red'); % Glider -5


% Set labels and title with LaTeX formatting
xlabel('$C_D$ (Induced Drag Coefficient)', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('$C_L$ (Lift Coefficient)', 'Interpreter', 'latex', 'FontSize', 14);
title('Lift-Drag Polar Curve', 'Interpreter', 'latex', 'FontSize', 16);

% Adjust axis for better visualization
axis tight;

legend('Cessna wing', 'Glider wing','Cessna with tail 0$^\circ$', 'Glider with tail 0$^\circ$', ...
    'Cessna with tail -2$^\circ$', 'Glider with tail -2$^\circ$', 'Interpreter', 'latex', 'FontSize', 12, 'Location', 'best');

exportgraphics(gcf, 'polar.pdf', 'ContentType', 'vector');

 %% Cl Alpha
figure();
plot(rad2deg(alpha_range), cL(:,1), '-bo', 'LineWidth', 1.5, 'MarkerSize', 6); % cessna
clAlphaC = polyfit(alpha_range, cL(:,1), 1);
hold on
plot(rad2deg(alpha_range), cL(:,2), '-ro', 'LineWidth', 1.5, 'MarkerSize', 6); % glider
plot(rad2deg(alpha_range), alpha_range*2*pi, 'g', 'LineWidth', 0.5);
grid on
plot(rad2deg(alpha_range), cL(:,3), '--b^', 'LineWidth', 1.5, 'MarkerSize', 6, 'MarkerFaceColor', 'blue'); % Cessna 0
plot(rad2deg(alpha_range), cL(:,4), '--r^', 'LineWidth', 1.5, 'MarkerSize', 6, 'MarkerFaceColor', 'red'); % Glider 0
plot(rad2deg(alpha_range), cL(:,5), '--g^', 'LineWidth', 1.5, 'MarkerSize', 6, 'MarkerFaceColor', 'blue'); % Cessna -5
plot(rad2deg(alpha_range), cL(:,6), '--m^', 'LineWidth', 1.5, 'MarkerSize', 6, 'MarkerFaceColor', 'red'); % Glider -5



% Set labels and title with LaTeX formatting
xlabel('$Alpha$ (AoA)', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('$C_L$ (Lift Coefficient)', 'Interpreter', 'latex', 'FontSize', 14);
title('Lift-AoA Polar Curve', 'Interpreter', 'latex', 'FontSize', 16);

% Adjust axis for better visualization
axis tight;

% Add a legend
legend('Cessna wing', 'Glider wing', 'Two Dimensional Thin Airfoil', ...
    'Cessna with tail 0$^\circ$', 'Glider with tail 0$^\circ$','Cessna with tail -2$^\circ$', ...
    'Glider with tail -2$^\circ$', 'Interpreter', 'latex', 'FontSize', 12, 'Location', 'best');


exportgraphics(gcf, 'cl_alpha.pdf', 'ContentType', 'vector');

%% Circulation distribution, AoA = 5°

alpha = deg2rad(5.0);
uInf = [cos(alpha), -sin(alpha), 0; sin(alpha), cos(alpha), 0; 0, 0, 1] * [1.0; 0.0; 0.0];
wingC = buildLinearSystem(uInf, wingC);
wingC = postprocessing(wingC);
wingG = buildLinearSystem(uInf, wingG);
wingG = postprocessing(wingG);


S = zeros(wingC.discretize(2));
thetaC = linspace(acos(- 2 * wingC.s(1)), acos(- 2 * wingC.s(end)), wingC.discretize(2));
for n = 1:wingC.discretize(2)
    S(:, n) = sin(n * thetaC);
end
a = S \ wingC.gammaDistribution'; % Fit Coefficients
% Build series
gammaEllipticC = S(:,1) * a(1);

S = zeros(wingG.discretize(2));
thetaG = linspace(acos(- 2 * wingG.s(1)), acos(- 2 * wingG.s(end)), wingG.discretize(2));
for n = 1:wingG.discretize(2)
    S(:, n) = sin(n * thetaG);
end
a = S \ wingG.gammaDistribution'; % Fit Coefficients
% Build series
gammaEllipticG = S(:,1) * a(1);

figure();
plot(wingC.s, wingC.gammaDistribution, '-bo', 'LineWidth', 1.5, 'MarkerSize', 6); % cessna
grid on;
hold on
plot(wingG.s, wingG.gammaDistribution, '-ro', 'LineWidth', 1.5, 'MarkerSize', 6); % glider
plot(-0.5 * cos(thetaC), gammaEllipticC, '-b', 'LineWidth', 0.5, 'MarkerSize', 6);
plot(-0.5 * cos(thetaG), gammaEllipticG, '-r', 'LineWidth', 0.5, 'MarkerSize', 6);
scatter([-0.5; 0.5], [0; 0]);


% Set labels and title with LaTeX formatting
xlabel('$s$ (Span Fraction)', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('$\Gamma$ (Circulation)', 'Interpreter', 'latex', 'FontSize', 14);
title('Circulation distribution', 'Interpreter', 'latex', 'FontSize', 16);

% Adjust axis for better visualization
axis tight;

% Optional: Add a legend if needed
legend('Cessna','Glider', 'Cessna Elliptic distribution','Glider Elliptic distribution',  'Interpreter', 'latex', 'FontSize', 12, 'Location', 'best');

set(gcf, 'Units', 'centimeters', 'Position', [0, 0, 20, 8]);
exportgraphics(gcf, 'circulation.pdf', 'ContentType', 'vector');


