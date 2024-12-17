clear; clc; close all
% Analysis of Cessna 172 Skyhawk wing, for method specifications read mainK.m

addpath('../');


%% Define element parametrically

% data
wing.ID = 1;
wing.xOffset = [0.0; 0.0; 0.0];
wing.rootChord = 1.68;
wing.span = 11.0;
wing.dihedral = deg2rad(3.0);
wing.sweep = deg2rad(0.0);
wing.taper = 0.6369;
wing.firstTaper = 0.418;
wing.twistPrime = deg2rad(2);
wing.airfoilCoefficients = [17.3648404437696	-87.1283457681812	185.023337210179	-216.177818020936	151.422278341809	-65.2483351673165	17.2902703811582	-2.77635906475406	0.150449092907973	0.0781532987533795	0.00158968881305051];


% extra
wing.tipChord = wing.rootChord * wing.taper;
wing.MGC = (wing.rootChord + wing.tipChord) / 2.0 * (1 - wing.firstTaper) + wing.rootChord * wing.firstTaper; % Mean Geometric Chord
wing.S = wing.MGC * wing.span;
wing.AR = wing.span^2 / wing.S;

% define precision
wing.discretize = [15; 50]; % singularities in [chord direction; spanwise direction]


% data
tail.ID = 2;
tail.xOffset = [4.88; 0.0; 0.0];
tail.rootChord = 0.94;
tail.span = 3.45;
tail.dihedral = deg2rad(0.0);
tail.sweep = deg2rad(0.0);
tail.taper = 0.4894;
tail.firstTaper = 0.6;
tail.twistPrime = deg2rad(0.0);
tail.airfoilCoefficients = [0.0; 0.0; 0.0];


% extra
tail.tipChord = tail.rootChord * tail.taper;
tail.MGC = (tail.rootChord + tail.tipChord) / 2.0 * (1 - tail.firstTaper) + tail.rootChord * tail.firstTaper; % Mean Geometric Chord
tail.S = tail.MGC * tail.span;
tail.AR = tail.span^2 / tail.S;

% define precision
tail.discretize = [5; 10]; % singularities in [chord direction; spanwise direction]


%% Build Elements
plotFlag = true;
wing = buildElement(wing, plotFlag);
tail = buildElement(tail, plotFlag);

%% Compute 
alpha_range = deg2rad(-5:2:12);
cL = zeros(size(alpha_range));
cD = cL;
for i = 1:length(alpha_range)
    alpha = alpha_range(i);
    uInf = [cos(alpha), -sin(alpha), 0; sin(alpha), cos(alpha), 0; 0, 0, 1] * [1.0; 0.0; 0.0];
    wing = buildLinearSystem(uInf, wing, tail);
    wing = postprocessing(wing);
    cL(i) = wing.cL;
    cD(i) = wing.cD;
end

alpha = deg2rad(2.0);
uInf = [cos(alpha), -sin(alpha), 0; sin(alpha), cos(alpha), 0; 0, 0, 1] * [1.0; 0.0; 0.0];
wing = buildLinearSystem(uInf, wing, tail);
wing = postprocessing(wing);
gammaD = wing.gammaDistribution';
s = wing.s;

%% Polar
figure(6);
plot(cD, cL, '-o', 'LineWidth', 1.5, 'MarkerSize', 6);
grid on;

% Set labels and title with LaTeX formatting
xlabel('$C_D$ (Induced Drag Coefficient)', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('$C_L$ (Lift Coefficient)', 'Interpreter', 'latex', 'FontSize', 14);
title('Lift-Drag Polar Curve', 'Interpreter', 'latex', 'FontSize', 16);

% Adjust axis for better visualization
axis tight;

% Optional: Add a legend if needed
legend('Simulation Data', 'Interpreter', 'latex', 'FontSize', 12, 'Location', 'best');

saveas(gcf, 'cessna_polar_wing.pdf');

%% Cl Alpha
figure(7);
plot(rad2deg(alpha_range), cL, '-o', 'LineWidth', 1.5, 'MarkerSize', 6);
grid on;

% Set labels and title with LaTeX formatting
xlabel('$C_D$ (AoA)', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('$C_L$ (Lift Coefficient)', 'Interpreter', 'latex', 'FontSize', 14);
title('Lift-AoA Polar Curve', 'Interpreter', 'latex', 'FontSize', 16);

% Adjust axis for better visualization
axis tight;

% Optional: Add a legend if needed
legend('Simulation Data', 'Interpreter', 'latex', 'FontSize', 12, 'Location', 'best');

saveas(gcf, 'cessna_alpha.pdf');

%% Circulation distribution, AoA = 2°

N = length(gammaD); % Sine functions number
S = zeros(length(s), N);
for n = 1:N
    S(:, n) = sin(n * pi * (s + 0.5));
end

a = S \ gammaD; % Fit Coefficients

% Build series
gammaFit = S * a;
gammaElliptic = S(:,1) * a(1);

figure(8);
plot(s, gammaD, '-o', 'LineWidth', 1.5, 'MarkerSize', 6);
grid on;
hold on
plot(s, gammaFit, 'LineWidth', 1.5, 'MarkerSize', 6);
plot(s, gammaElliptic, 'LineWidth', 1.5, 'MarkerSize', 6);
scatter([-0.5; 0.5], [0; 0]);


% Set labels and title with LaTeX formatting
xlabel('$s$ (Span Fraction)', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('$\Gamma$ (Circulation)', 'Interpreter', 'latex', 'FontSize', 14);
title('Circulation distribution', 'Interpreter', 'latex', 'FontSize', 16);

% Adjust axis for better visualization
axis tight;

% Optional: Add a legend if needed
legend('Simulation Data','Sine fitting', 'Elliptic distribution',  'Interpreter', 'latex', 'FontSize', 12, 'Location', 'best');

saveas(gcf, 'cessna_gamma_distribution.pdf');

