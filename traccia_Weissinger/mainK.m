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
wing.twistPrime = deg2rad(20);
% airfoils are described by a fitting of the mean line over three sine
% functions y = A*sin(pi *x) + B*sin(2pi *x) + C*sin(3pi *x), this can be
% easily changed with real airfoil mean line points,
wing.airfoilCoefficients = [0.05; 0.03; 0.01];


% extra
wing.tipChord = wing.rootChord * wing.taper;
wing.MGC = (wing.rootChord + wing.tipChord) / 2.0; % Mean Geometric Chord
wing.S = wing.MGC * wing.span;
wing.AR = wing.span^2 / wing.S;

% define precision
wing.discretize = [20; 50]; % singularities in [chord direction; spanwise direction]


tic
wing = buildElement(wing);
t = toc


