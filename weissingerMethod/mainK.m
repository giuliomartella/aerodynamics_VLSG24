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
wing.discretize = [15; 50]; % singularities in [chord direction; spanwise direction]


% data
tail.ID = 2;
tail.xOffset = [10.0; 0.0; 0.0];
tail.rootChord = 0.50;
tail.span = 3.0;
tail.dihedral = deg2rad(0.0);
tail.sweep = deg2rad(5.0);
tail.taper = 0.5;
tail.twistPrime = deg2rad(0);
tail.airfoilCoefficients = [0.0; 0.0; 0.0];


% extra
tail.tipChord = tail.rootChord * tail.taper;
tail.MGC = (tail.rootChord + tail.tipChord) / 2.0; % Mean Geometric Chord
tail.S = tail.MGC * tail.span;
tail.AR = tail.span^2 / tail.S;

% define precision
tail.discretize = [5; 10]; % singularities in [chord direction; spanwise direction]


%% Build Elements
plotFlag = true;
tic
wing = buildElement(wing, plotFlag);
tail = buildElement(tail, plotFlag);

elapsedTime0 = toc;
disp(['Elapsed time for construction: ', num2str(elapsedTime0), ' seconds.']);


%% Define far field conditions
uInf = [1.0; 0.0; 0.0];
alpha = deg2rad(1.0);
uInf = [cos(alpha), -sin(alpha), 0; sin(alpha), cos(alpha), 0; 0, 0, 1] * uInf;



%% Build Linear System
% 
% [wing, tail] = buildLinearSystem(uInf, wing, tail);
% elapsedTime1 = toc;
% disp(['Elapsed time for solving linear system: ', num2str(elapsedTime1 - elapsedTime0), ' seconds.']);

%% Circulation, Lift, Drag
% 
% wing = postprocessing(wing);
% tail = postprocessing(tail);



%% Loop for Polar Plot

elapsedTime2 = toc;

alpha_range = deg2rad(-5:2:12);  
polarPlot(alpha_range, wing, tail);
% polarPlot(alpha_range, wing);

elapsedTime3 = toc;
disp(['Elapsed time for computing polars: ', num2str(elapsedTime3 - elapsedTime2), ' seconds.']);
disp(['Average elapsed time for polar point: ', num2str((elapsedTime3 - elapsedTime2)/ length(alpha_range)) , ' seconds.']);
