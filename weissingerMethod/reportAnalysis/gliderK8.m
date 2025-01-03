function [wingG, tailG] = gliderK8()
%% Define element parametrically

% data
wingG.ID = 1;
wingG.xOffset = [0.0; 0.0; 0.0];
wingG.rootChord = 1.05;
wingG.span = 15.0;
wingG.dihedral = deg2rad(2.0);
wingG.sweep = deg2rad(0.0);
wingG.taper = 0.3619;
wingG.firstTaper = 0.415;
wingG.twistPrime = deg2rad(0);
wingG.twistZero = deg2rad(0.0);
wingG.airfoilCoefficients = - [-109.381022564051	548.763980581295	-1174.91212765915	1402.64294899680	-1024.59095207000	473.802207095480	-139.597200898110	26.3778326749395	-3.57684985871406	0.469337303906058	0.00147990272037135];


% extra
wingG.tipChord = wingG.rootChord * wingG.taper;
wingG.MGC = (wingG.rootChord + wingG.tipChord) / 2.0 * (1 - wingG.firstTaper) + wingG.rootChord * wingG.firstTaper; % Mean Geometric Chord
wingG.S = wingG.MGC * wingG.span;
wingG.AR = wingG.span^2 / wingG.S;

% define precision
wingG.discretize = [30; 30]; % singularities in [chord direction; spanwise direction]



% data
tailG.ID = 2;
tailG.xOffset = [5.60; -0.5; 0.0];
tailG.rootChord = 0.65;
tailG.span = 3.20;
tailG.dihedral = deg2rad(0.0);
tailG.sweep = deg2rad(0.0);
tailG.taper = 0.4923;
tailG.firstTaper = 0.95;
tailG.twistPrime = deg2rad(0.0);
tailG.twistZero = deg2rad(0.0);
tailG.airfoilCoefficients = [0.0; 0.0; 0.0];


% extra
tailG.tipChord = tailG.rootChord * tailG.taper;
tailG.MGC = (tailG.rootChord + tailG.tipChord) / 2.0 * (1 - tailG.firstTaper) + tailG.rootChord * tailG.firstTaper; % Mean Geometric Chord
tailG.S = tailG.MGC * tailG.span;
tailG.AR = tailG.span^2 / tailG.S;

% define precision
tailG.discretize = [20; 30]; % singularities in [chord direction; spanwise direction]



end