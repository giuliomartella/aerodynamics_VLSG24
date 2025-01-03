function wing = buildElement(wing, plotFlag)
% Each aerodynamic surface is described as a vortex distribution. 
% Airfoils are represented by a n-order polynomial.

% singularities in [chord direction; spanwise direction]



%% Control Points
wing.controlPoint = zeros(wing.discretize(1), wing.discretize(2), 3);
s = linspace(-0.5, 0.5, wing.discretize(2)+1);
s = s(1: wing.discretize(2)) + 1 / wing.discretize(2) / 2; % abscissa
wing.chordDistribution = (1 - (abs(s)/0.5 - wing.firstTaper) * wing.taper) * wing.rootChord;
wing.chordDistribution(wing.chordDistribution > wing.rootChord) = wing.rootChord;
wing.twistDistribution = wing.twistPrime * abs(s)/0.5 - wing.twistZero;

wing.quarter = zeros(wing.discretize(2), 3);
wing.quarter(:, 1) = wing.chordDistribution * 0.25;
wing.quarter(:, 3) = s * wing.span * cos(wing.dihedral) * cos(wing.sweep);

for i = 1:wing.discretize(2)
    % define planar distribution 
    wing.controlPoint(:, i, 1) = linspace(wing.chordDistribution(i) / wing.discretize(1) * 0.75, wing.chordDistribution(i) * (1 - 1 / wing.discretize(1) * 0.25), wing.discretize(1));
    wing.controlPoint(:, i, 2) = findY(wing.airfoilCoefficients, wing.chordDistribution(i), wing.twistDistribution(i), wing.controlPoint(:, i, 1));
    wing.controlPoint(:, i, 3) = s(i) * wing.span * cos(wing.dihedral) * cos(wing.sweep);

    wing.quarter(i, 2) = findY(wing.airfoilCoefficients, wing.chordDistribution(i), wing.twistDistribution(i), wing.quarter(i, 1));
end

% add third dimension to x and y
wing.controlPoint(:, :, 1) = wing.controlPoint(:, :, 1) + abs(wing.controlPoint(:, :, 3)) * sin(wing.sweep);
wing.controlPoint(:, :, 2) = wing.controlPoint(:, :, 2) + abs(wing.controlPoint(:, :, 3)) * sin(wing.dihedral);
wing.quarter(:, 1) =  wing.quarter(:, 1) + abs(wing.quarter(:, 3)) * sin(wing.sweep);
wing.quarter(:, 2) =  wing.quarter(:, 2) + abs(wing.quarter(:, 3)) * sin(wing.dihedral);

% add offset
wing.controlPoint = wing.controlPoint + reshape(wing.xOffset, 1, 1, 3);
wing.quarter = wing.quarter + wing.xOffset';

wing.s = s;

%% Normal vector 


wing.normal = zeros(wing.discretize(1), wing.discretize(2), 3);

% We may now take advantage of the mean line analytic description.
yPrime = zeros(wing.discretize(1), wing.discretize(2));
quarterYPrime = zeros(1, wing.discretize(2));
for i = 1:wing.discretize(2)
    yPrime(:, i) = findYPrime(wing.airfoilCoefficients, wing.chordDistribution(i), wing.twistDistribution(i), wing.controlPoint(:, i, 1));
    quarterYPrime(i) = findYPrime(wing.airfoilCoefficients, wing.chordDistribution(i), wing.twistDistribution(i), wing.quarter(i, 1));
end
% Before considering dihedral effect the tangent versor is tau =
% [1; yPrime], so the normal is n = [-yPrime; 1]
wing.normal(:,:, 1) = - yPrime ./ sqrt(yPrime.^2 +1) * cos(wing.dihedral);
wing.normal(:,:, 2) = 1 ./ sqrt(yPrime.^2 +1) * cos(wing.dihedral);

% Now consider dihedral
wing.normal(:, :, 3) = sin(wing.dihedral) .* repmat(sign(-s(:))', size(wing.normal, 1), 1);

% Normal to quarter chord points
wing.quarterNormal(:, 1) = - quarterYPrime ./ sqrt(quarterYPrime.^2 +1) * cos(wing.dihedral);
wing.quarterNormal(:, 2) = 1 ./ sqrt(quarterYPrime.^2 +1) * cos(wing.dihedral);
wing.quarterNormal(:, 3) = sin(wing.dihedral) .* sign(-s(:))';

%% Vortices defining vertex positions
% VFL: vertex forward left
% VteL : vertex on trailing edge left
% VBR: vertex backward right
% For simplicity and cost effctiveness vortex position is shaped on
% linearized mean line
spanwiseDelta = wing.span / wing.discretize(2) / 2.0; % scalar
addSpanDelta = zeros(wing.discretize(1), wing.discretize(2), 3);
addSpanDelta(:, :, 3) = ones(wing.discretize') * spanwiseDelta;

chordwiseDelta = wing.chordDistribution / wing.discretize(1) * 0.25; % vector
addChordDelta = zeros(wing.discretize(1), wing.discretize(2), 3);
addChordDelta(:, :, 1) = repmat(chordwiseDelta, wing.discretize(1), 1);

wing.VFL = wing.controlPoint - addChordDelta - addSpanDelta;
wing.VFR = wing.controlPoint - addChordDelta + addSpanDelta;

% Building the linear system the scrpit will add a distance in uInf
% direction.

wing.VBL = wing.controlPoint + addChordDelta - addSpanDelta;
wing.VBR = wing.controlPoint + addChordDelta + addSpanDelta;

wing.VteL = repmat(wing.VBL(end, :, :), wing.discretize(1), 1, 1) + addChordDelta;
wing.VteR = repmat(wing.VBR(end, :, :), wing.discretize(1), 1, 1) + addChordDelta;


%% functions
    function y = findY(airfoilCoefficients, chord, twist, cp)
        sD = cp ./ chord; % s dummy
        y = polyval(airfoilCoefficients, sD);
        twisted = [cos(twist), -sin(twist); sin(twist), cos(twist)] * [sD, y]';
        y = twisted(2, :);

    end

    function yPrime = findYPrime(airfoilCoefficients, chord, twist, cp)
        sD = cp ./ chord; % s dummy
        yPrime =  polyval(polyder(airfoilCoefficients), sD);
        twisted = [cos(twist), -sin(twist); sin(twist), cos(twist)] * [sD, yPrime]';
        yPrime = twisted(2, :);
    end

%% plot functions
% fancy plots to check wing shape, plot functions are AI generated
if plotFlag
    plotWingNormals(wing)
    scatterControlPoints(wing)
    plotWingSurface(wing)
    plotWingWithVortices(wing)
    plotQuarterChordWithNormals(wing)
end

    function scatterControlPoints(wing)
        % SCATTERCONTROLPOINTS Visualizes the 3D control points of a wing
        %
        % This function takes the wing structure containing the control points
        % and plots them in 3D space using a scatter plot.
        %
        % Inputs:
        %   - wing: A structure with a field `controlPoint`, which is a 3D array
        %           of dimensions [N_chord x N_span x 3], representing the
        %           coordinates of control points in 3D space.

    % Extract control points
    % `controlPoint` is a 3D array: [x, y, z] for each control point
    controlPoints = wing.controlPoint;

    % Reshape the control points for easier plotting
    % The 3rd dimension contains [x, y, z] coordinates
    x = controlPoints(:, :, 1);
    y = controlPoints(:, :, 2);
    z = controlPoints(:, :, 3);

    % Flatten the matrices to vectors for scatter plotting
    x = x(:); % x-coordinates
    y = y(:); % y-coordinates
    z = z(:); % z-coordinates

    % Scatter the control points in 3D space
    figure(1); % Open a new figure
    scatter3(x, y, z, 20, 'filled'); % Scatter plot, marker size = 50
    grid on; % Enable grid for better visualization
    axis equal; % Use equal scaling for all axes
    xlabel('X (Chord direction)'); % Label for x-axis
    ylabel('Y (Vertical direction)'); % Label for y-axis
    zlabel('Z (Spanwise direction)'); % Label for z-axis
    title('3D Scatter Plot of Control Points'); % Plot title

    % Add markers or lines to indicate the wing's structure if needed
    % Uncomment the following line to add a connection between points
    % hold on; plot3(x, y, z, 'k-', 'LineWidth', 0.5); % Example connection
    hold on

end

function plotWingSurface(wing)
% PLOTWINGSURFACE Visualizes the wing surface in 3D using a surface plot
%
% This function creates a 3D representation of the wing by plotting
% the control points as a continuous surface.
%
% Inputs:
%   - wing: A structure with a field `controlPoint`, which is a 3D array 
%           of dimensions [N_chord x N_span x 3], representing the 
%           coordinates of control points in 3D space.

    % Extract control points
    % `controlPoint` is a 3D array: [x, y, z] for each control point
    controlPoints = wing.controlPoint;

    % Extract the individual coordinate matrices
    x = controlPoints(:, :, 1); % Chordwise direction
    y = controlPoints(:, :, 2); % Vertical displacement (airfoil shape)
    z = controlPoints(:, :, 3); % Spanwise direction

    % Plot the surface
    figure(2); % Open a new figure
    surf(z, x, y, 'FaceColor', 'interp', 'EdgeColor', 'none'); 
    % `surf` assumes X, Y, Z -> we map Span (z), Chord (x), Vertical (y)

    % Aesthetic improvements
    colormap jet; % Use a color map for surface shading
    colorbar; % Add a color bar to visualize the height variation
    grid on; % Enable grid
    axis equal; % Use equal scaling for all axes
    xlabel('Spanwise Direction (Z)'); % Label for x-axis
    ylabel('Chord Direction (X)'); % Label for y-axis
    zlabel('Vertical Displacement (Y)'); % Label for z-axis
    title('3D Wing Surface Visualization'); % Plot title
    view(3); % Set 3D view angle
    hold on
end

function plotWingNormals(wing)
% PLOTWINGNORMALS Visualizes the normal vectors of the wing surface
%
% This function plots the wing surface along with its normal vectors 
% using quiver3 for clear representation.
%
% Inputs:
%   - wing: A structure with fields:
%       * controlPoint: 3D array of control points [N_chord x N_span x 3]
%       * normal: 3D array of normal vectors [N_chord x N_span x 3]

    % Extract control points
    controlPoints = wing.controlPoint;
    x = controlPoints(:, :, 1); % Chordwise direction
    y = controlPoints(:, :, 2); % Vertical displacement (airfoil shape)
    z = controlPoints(:, :, 3); % Spanwise direction

    % Extract normal vectors
    normals = wing.normal;
    nx = normals(:, :, 1); % Normal vector X component
    ny = normals(:, :, 2); % Normal vector Y component
    nz = normals(:, :, 3); % Normal vector Z component

    % Plot the wing surface
    figure(3);
    surf(z, x, y, 'FaceColor', 'interp', 'EdgeColor', 'none'); 
    hold on;

    % Add the normal vectors
    % Flatten the arrays for quiver3
    x_flat = x(:);
    y_flat = y(:);
    z_flat = z(:);
    nx_flat = nx(:);
    ny_flat = ny(:);
    nz_flat = nz(:);

    % Scale factor for normal vector length
    scaleFactor = 0.1 * max(wing.span, max(wing.chordDistribution(:))); 

    % Plot normal vectors as arrows
    quiver3(z_flat, x_flat, y_flat, ...
            scaleFactor * nz_flat, scaleFactor * nx_flat, scaleFactor * ny_flat, ...
            'k', 'LineWidth', 1, 'MaxHeadSize', 0.5);

    % Aesthetic improvements
    colormap jet; % Use a color map for surface shading
    colorbar; % Add a color bar to visualize the height variation
    grid on; % Enable grid
    axis equal; % Use equal scaling for all axes
    xlabel('Spanwise Direction (Z)'); % Label for x-axis
    ylabel('Chord Direction (X)'); % Label for y-axis
    zlabel('Vertical Displacement (Y)'); % Label for z-axis
    title('3D Wing Surface with Normal Vectors'); % Plot title
    view(3); % Set 3D view angle

end

function plotWingWithVortices(wing)
% PLOTWINGWITHVORTICES Visualizes the wing surface and vortex vertices in 3D
%
% This function creates a 3D representation of the wing by plotting
% the control points as a continuous surface and overlays the vortex 
% vertices (VFL, VFR, VBL, VBR) as scatter points.
%
% Inputs:
%   - wing: A structure with fields:
%       * controlPoint: 3D array of control points [N_chord x N_span x 3]
%       * VFL, VFR, VBL, VBR: 3D arrays of vortex vertices 
%                             [N_chord x N_span x 3]

    % Extract control points
    controlPoints = wing.controlPoint;
    x = controlPoints(:, :, 1); % Chordwise direction
    y = controlPoints(:, :, 2); % Vertical displacement (airfoil shape)
    z = controlPoints(:, :, 3); % Spanwise direction

    % Plot the wing surface
    figure(4);
    surf(z, x, y, 'FaceColor', 'interp', 'EdgeColor', 'none'); 
    hold on;

    % Overlay vortex vertices
    scatterVortexVertices(wing.VFL, 'r', 'VFL');
    scatterVortexVertices(wing.VFR, 'g', 'VFR');
    % scatterVortexVertices(wing.VBL, 'b', 'VBL');
    % scatterVortexVertices(wing.VBR, 'k', 'VBR');

    % Aesthetic improvements
    colormap jet; % Use a color map for surface shading
    colorbar; % Add a color bar to visualize height variation
    grid on; % Enable grid
    axis equal; % Use equal scaling for all axes
    xlabel('Spanwise Direction (Z)'); % Label for x-axis
    ylabel('Chord Direction (X)'); % Label for y-axis
    zlabel('Vertical Displacement (Y)'); % Label for z-axis
    title('3D Wing Surface with Vortex Vertices'); % Plot title
    view(3); % Set 3D view angle


    % Nested function to scatter vortex vertices
    function scatterVortexVertices(vertices, color, label)
        % SCATTERVORTEXVERTICES Plots vortex vertices as scatter points
        %
        % Inputs:
        %   - vertices: 3D array of vertex coordinates [N_chord x N_span x 3]
        %   - color: Color for the scatter points (e.g., 'r' for red)
        %   - label: String label for legend (optional)

        % Reshape vertices for scatter plot
        vx = vertices(:, :, 1); % X-coordinates
        vy = vertices(:, :, 2); % Y-coordinates
        vz = vertices(:, :, 3); % Z-coordinates

        scatter3(vz(:), vx(:), vy(:), 20, color, 'filled'); % Scatter points
    end
end


function plotQuarterChordWithNormals(wing)
% PLOTQUARTERCHORDWITHNORMALS Visualizes the quarter chord points and their normals in 3D.
%
% This function plots the quarter chord points on the wing surface and overlays
% their normals for visualization, ensuring consistency with the wing surface plot.
%
% Inputs:
%   - wing: A structure containing wing geometry data, including:
%       * wing.controlPoint: 3D array [N_chord x N_span x 3] of surface points.
%       * wing.quarter: 3D array [N_chord x N_span x 3] of quarter chord points.
%       * wing.quarterNormal: 3D array [N_chord x N_span x 3] of normal vectors at quarter chord points.

    % Extract quarter chord points and their normals
    xQuarterChord = wing.quarter(:, 1); % X-coordinates
    yQuarterChord = wing.quarter(:, 2); % Y-coordinates
    zQuarterChord = wing.quarter(:, 3); % Z-coordinates

    normals = wing.quarterNormal; % Normal vectors at quarter chord points
    nx = normals(:, 1); % X-component of normals
    ny = normals(:, 2); % Y-component of normals
    nz = normals(:, 3); % Z-component of normals

    % Plot the wing surface (consistent with figure 2)
    figure(2);
    hold on;

    % Overlay quarter chord points
    plot3(zQuarterChord, xQuarterChord, yQuarterChord, 'o-', 'DisplayName', 'Quarter Chord Line');

    % Add normal vectors to the plot
    quiver3(zQuarterChord, xQuarterChord, yQuarterChord, ...
            nz, nx, ny, ...
            0.5, 'r', 'LineWidth', 1.5, 'DisplayName', 'Normals');

    % Aesthetic improvements
    xlabel('Spanwise Direction (Z)');
    ylabel('Chord Direction (X)');
    zlabel('Vertical Displacement (Y)');
    grid on;
    axis equal;
    legend('show');
    title('Quarter Chord Points and Normals');
    hold off;
end


end