function [wingC, tailC] = cessna()
%% Define element parametrically

% data
wingC.ID = 1;
wingC.xOffset = [0.0; 0.0; 0.0];
wingC.rootChord = 1.68;
wingC.span = 11.0;
wingC.dihedral = deg2rad(3.0);
wingC.sweep = deg2rad(0.0);
wingC.taper = 0.6369;
wingC.firstTaper = 0.418;
wingC.twistPrime = deg2rad(2);
wingC.airfoilCoefficients = [17.3648404437696	-87.1283457681812	185.023337210179	-216.177818020936	151.422278341809	-65.2483351673165	17.2902703811582	-2.77635906475406	0.150449092907973	0.0781532987533795	0.00158968881305051];


% extra
wingC.tipChord = wingC.rootChord * wingC.taper;
wingC.MGC = (wingC.rootChord + wingC.tipChord) / 2.0 * (1 - wingC.firstTaper) + wingC.rootChord * wingC.firstTaper; % Mean Geometric Chord
wingC.S = wingC.MGC * wingC.span;
wingC.AR = wingC.span^2 / wingC.S;

% define precision
wingC.discretize = [21; 41]; % singularities in [chord direction; spanwise direction]


% data
tailC.ID = 2;
tailC.xOffset = [4.88; 0.0; 0.0];
tailC.rootChord = 0.94;
tailC.span = 3.45;
tailC.dihedral = deg2rad(0.0);
tailC.sweep = deg2rad(0.0);
tailC.taper = 0.4894;
tailC.firstTaper = 0.6;
tailC.twistPrime = deg2rad(0.0);
tailC.airfoilCoefficients = [0.0; 0.0; 0.0];


% extra
tailC.tipChord = tailC.rootChord * tailC.taper;
tailC.MGC = (tailC.rootChord + tailC.tipChord) / 2.0 * (1 - tailC.firstTaper) + tailC.rootChord * tailC.firstTaper; % Mean Geometric Chord
tailC.S = tailC.MGC * tailC.span;
tailC.AR = tailC.span^2 / tailC.S;

% define precision
tailC.discretize = [6; 11]; % singularities in [chord direction; spanwise direction]


end