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
wing.twistPrime = deg2rad(10);
% airfoils are described by a fitting of the mean line over three sine
% functions y = A*sin(pi *x) + B*sin(2pi *x) + C*sin(3pi *x), this can be
% easily changed with real airfoil mean line points,
wing.airfoilCoefficients = [0.007654514121448;-0.002202436502870;9.065825636278312e-04];


% extra
wing.tipChord = wing.rootChord * wing.taper;
wing.MGC = (wing.rootChord + wing.tipChord) / 2.0; % Mean Geometric Chord
wing.S = wing.MGC * wing.span;
wing.AR = wing.span^2 / wing.S;

% define precision
wing.discretize = [10; 30]; % singularities in [chord direction; spanwise direction]


% data
tail.ID = 2;
tail.xOffset = [5.0; 0.0; 0.0];
tail.rootChord = 0.50;
tail.span = 3.0;
tail.dihedral = deg2rad(0.0);
tail.sweep = deg2rad(5.0);
tail.taper = 0.5;
tail.twistPrime = deg2rad(2);
% airfoils are described by a fitting of the mean line over three sine
% functions y = A*sin(pi *x) + B*sin(2pi *x) + C*sin(3pi *x), this can be
% easily changed with real airfoil mean line points,
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
alpha = deg2rad(2.0);
uInf = [cos(alpha), -sin(alpha), 0; sin(alpha), cos(alpha), 0; 0, 0, 1] * uInf;



%% Build Linear System

gamma = buildLinearSystem(uInf, wing, tail);
elapsedTime1 = toc;
disp(['Elapsed time for solving linear system: ', num2str(elapsedTime1 - elapsedTime0), ' seconds.']);






