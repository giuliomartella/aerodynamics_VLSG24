clear; clc; close all
% K stands for KISS: "keep it simple and straight forward"


%% Define element parametrically
% Sweep and Dihedral angles are spanwise constant and use leading edge as
% reference. Taper ratio is spanwise constant and defined as TR =
% tipChord/rootChord. Twist is described as linear, the coefficient
% represent derivative over adimensional span abscissa.

% data
wing.ID = 1;
wing.xOffset = [0.0; 0.0; 0.0];
wing.rootChord = 1.0;
wing.span = 10.0;
wing.dihedral = deg2rad(2.0);
wing.sweep = deg2rad(5.0);
wing.taper = 0.5;
wing.twistPrime = deg2rad(3);
% airfoils are described by a fitting of the mean line over a second order
% polynomial
wing.airfoilCoefficients = [0.0; 0.0; 0.0];


% extra
wing.tipChord = wing.rootChord * wing.taper;
wing.MGC = (wing.rootChord + wing.tipChord) / 2.0; % Mean Geometric Chord
wing.S = wing.MGC * wing.span;
wing.AR = wing.span^2 / wing.S;

% define precision
wing.discretize = [10; 40]; % singularities in [chord direction; spanwise direction]


% data
tail.ID = 2;
tail.xOffset = [5.0; 0.0; 0.0];
tail.rootChord = 0.50;
tail.span = 3.0;
tail.dihedral = deg2rad(2.0);
tail.sweep = deg2rad(5.0);
tail.taper = 0.5;
tail.twistPrime = deg2rad(2);
tail.airfoilCoefficients = [0.0; 0.0; 0.0];


% extra
tail.tipChord = tail.rootChord * tail.taper;
tail.MGC = (tail.rootChord + tail.tipChord) / 2.0; % Mean Geometric Chord
tail.S = tail.MGC * tail.span;
tail.AR = tail.span^2 / tail.S;

% define precision
tail.discretize = [5; 10]; % singularities in [chord direction; spanwise direction]


%% Build Elements
tic
wing = buildElement(wing);
tail = buildElement(tail);

elapsedTime0 = toc;
disp(['Elapsed time for construction: ', num2str(elapsedTime0), ' seconds.']);


%% Define far field conditions
uInf = [1.0; 0.0; 0.0];
alpha = deg2rad(1.0);
uInf = [cos(alpha), -sin(alpha), 0; sin(alpha), cos(alpha), 0; 0, 0, 1] * uInf;



%% Build Linear System

[wing, tail] = buildLinearSystem(uInf, wing, tail);
elapsedTime1 = toc;
disp(['Elapsed time for solving linear system: ', num2str(elapsedTime1 - elapsedTime0), ' seconds.']);

%% Circulation, Lift, Drag
% 
% wing = postprocessing(wing);
% tail = postprocessing(tail);



%% Loop for Polar Plot

elapsedTime2 = toc;

alpha_range = deg2rad(-5:2:14);  
cl_values_wing = zeros(size(alpha_range));
cl_values_tail = zeros(size(alpha_range));
cl2d_wing = zeros(length(alpha_range), wing.discretize(2));  
cl2d_tail = zeros(length(alpha_range), tail.discretize(2));  

% Define the span range for wing and tail
span_wing = linspace(-0.5 * wing.span, 0.5 * wing.span, wing.discretize(2));  
span_tail = linspace(-0.5 * tail.span, 0.5 * tail.span, tail.discretize(2));

for i = 1:length(alpha_range)
    alpha = alpha_range(i);
    uInf = [cos(alpha), -sin(alpha), 0; sin(alpha), cos(alpha), 0; 0, 0, 1] * [1.0; 0.0; 0.0];
    [wing, tail] = buildLinearSystem(uInf, wing, tail);
    wing = postprocessing(wing);
    tail = postprocessing(tail);
    
    cl_values_wing(i) = wing.cL;  
    cl_values_tail(i) = tail.cL;
    cl2d_wing(i, :) = wing.cL2D;
    cl2d_tail(i, :) = tail.cL2D;
end

figure;

% Polar plot for the wing
subplot(2, 2, 1);
plot(alpha_range * 180/pi, cl_values_wing, 'o-', 'LineWidth', 2);
xlabel('Angle of Attack (degrees)');
ylabel('Lift Coefficient (C_L)');
title('Wing Polar Plot');
grid on;

% Polar plot for the tail
subplot(2, 2, 2);
plot(alpha_range * 180/pi, cl_values_tail, 'o-', 'LineWidth', 2);
xlabel('Angle of Attack (degrees)');
ylabel('Lift Coefficient (C_L)');
title('Tail Polar Plot');
grid on;

% Lift distribution for the wing (span-wise)
subplot(2, 2, 3);
for i = 1:length(alpha_range)
    plot(span_wing, cl2d_wing(i, :), 'LineWidth', 2);
    hold on;
end
xlabel('Span (m)');
ylabel('Lift Distribution (C_L2D)');
title('Wing Lift Distribution');
grid on;

% Lift distribution for the tail (span-wise)
subplot(2, 2, 4);
for i = 1:length(alpha_range)
    plot(span_tail, cl2d_tail(i, :), 'LineWidth', 2);
    hold on;
end
xlabel('Span (m)');
ylabel('Lift Distribution (C_L2D)');
title('Tail Lift Distribution');
grid on;


elapsedTime3 = toc;
disp(['Elapsed time for computing polars: ', num2str(elapsedTime3 - elapsedTime2), ' seconds.']);

