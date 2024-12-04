function wing = buildElement(wing, discretize)
% Each aerodynamic surface is described as a vortex distribution. 
% Airfoils are represented by a sine series

% singularities in [chord direction; spanwise direction]



%% Control Points
wing.controlPoint = zeros(discretize(1), discretize(2), 3);
x = linspace(-0.5, 0.5, discretize(2)+1);
x = x(1: discretize(2)) + 1 / discretize(2) / 2; % abscissa
wing.chordDistribution = (1- abs(x)/0.5 * wing.taper) * wing.rootChord;
wing.twistDistribution = wing.twistPrime * abs(x)/0.5;

for i = 1:discretize(2)
    % define planar distribution 
    wing.controlPoint(:, i, 1) = linspace(0, wing.chordDistribution(i), discretize(1));
    wing.controlPoint(:, i, 2) = findY(wing.airfoilCoefficients, wing.chordDistribution(i), wing.twistDistribution(i), discretize);
    wing.controlPoint(:, i, 3) = x(i) * wing.span * cos(wing.dihedral) * cos(wing.sweep);


end

% add third dimension to x and y
wing.controlPoint(:, :, 1) = wing.controlPoint(:, :, 1) + abs(wing.controlPoint(:, :, 3)) * sin(wing.sweep);
wing.controlPoint(:, :, 2) = wing.controlPoint(:, :, 2) + abs(wing.controlPoint(:, :, 3)) * sin(wing.dihedral);

% fancy plots to check wing shape, plot functions are AI generated
% scatterControlPoints(wing)
 plotWingSurface(wing)

%% Normal vector 


wing.normal = zeros(discretize(1), discretize(2), 3);

% We may now take advantage of the sine series description of the mean
% line. x-derivative will be obtained analytically.
yPrime = zeros(discretize(1), discretize(2));
for i = 1:discretize(2)
    yPrime(:, i) = findYPrime(wing.airfoilCoefficients, wing.chordDistribution(i), wing.twistDistribution(i), discretize);
end
% Before considering dihedral effect the tangent versor is tau =
% [1; yPrime], so the normal is n = [-yPrime; 1]
wing.normal(:,:, 1) = - yPrime ./ sqrt(yPrime.^2 +1) * cos(wing.dihedral);
wing.normal(:,:, 2) = 1 ./ sqrt(yPrime.^2 +1) * cos(wing.dihedral);

% Now consider dihedral
wing.normal(:,:, 3) = sin(wing.dihedral) .* sign(- x);

plotWingNormals(wing)








%% functions
    function y = findY(airfoilCoefficients, chordDistribution, twist, discretize)
        s = linspace(0, chordDistribution, discretize(1));
        y = airfoilCoefficients(1) * sin(1 * pi * s) ...
            + airfoilCoefficients(2) * sin(2 * pi * s)  ...
            + airfoilCoefficients(3) * sin(3 * pi * s);
        twisted = [cos(twist), -sin(twist); sin(twist), cos(twist)] * [s; y];
        y = twisted(2, :);

    end

    function yPrime = findYPrime(airfoilCoefficients, chordDistribution, twist, discretize)
        s = linspace(0, chordDistribution, discretize(1));
        yPrime = airfoilCoefficients(1) * cos(1 * pi * s) * 1 * pi ...
            + airfoilCoefficients(2) * cos(2 * pi * s) * 2 * pi ...
            + airfoilCoefficients(3) * cos(3 * pi * s) * 3 * pi;
        twisted = [cos(twist), -sin(twist); sin(twist), cos(twist)] * [s; yPrime];
        yPrime = twisted(2, :);
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
    figure; % Open a new figure
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
    figure; % Open a new figure
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
    figure;
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

    hold off;
end

end