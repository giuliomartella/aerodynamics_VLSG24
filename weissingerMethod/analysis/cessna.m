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
wingC.airfoilCoefficients = -[6953095849.62522	-41705878717.0905	89751168614.7455	-62976216426.9082	-43436941475.1917	77420233625.6567	-16159611478.3620	-12513098258.7942	19420850389.8692	-17935551589.3764	-36529797516.7625	41867351285.1980	6876558740.08507	416901458.244638	-2754388980.54345	-5849976232.15455	-16905370326.1947	-9983232397.21804	21637030424.7906	29759422518.9410	-8253321713.41345	-32225031040.4622	7070571675.04805	-18444815310.1994	34588676029.8893	950540000.402188	7180624381.41682	-23690512730.7859	-22152430753.6558	15700440104.8325	62675825986.7256	-66413711709.2122	-1312158333.09356	11114628842.5214	36839654139.9453	-49549699654.2755	11395033429.3859	24803453554.2620	-30987574581.6271	19582986085.5621	-8246066455.35853	2496444946.79758	-558646216.794745	93034211.9387504	-11457741.3792978	1025361.54320094	-64688.8128179092	2741.29645066055	-72.1588578634163	1.10672018699355	-0.00418445449216117];
% extra
wingC.tipChord = wingC.rootChord * wingC.taper;
wingC.MGC = (wingC.rootChord + wingC.tipChord) / 2.0 * (1 - wingC.firstTaper) + wingC.rootChord * wingC.firstTaper; % Mean Geometric Chord
wingC.S = wingC.MGC * wingC.span;
wingC.AR = wingC.span^2 / wingC.S;

% define precision
wingC.discretize = [50; 40]; % singularities in [chord direction; spanwise direction]



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
tailC.discretize = [10; 20]; % singularities in [chord direction; spanwise direction]


end