function wing = buildElement(wing, discretize)
% Each aerodynamic surface is described as a vortex distribution. 
% Airfoils are represented by a sine series

% singularities in [chord direction; spanwise direction]



%% Control Points
nSingularities = discretize(1) * discretize(2);
wing.controlPoint = zeros(discretize(1), discretize(2), 3);
x = linspace(-0.5, 0.5, discretize(2)+1);
x = x(1: discretize(2)) + 1 / discretize(2) / 2; % abscissa
wing.chordDistribution = (1- abs(x)/0.5 * wing.taper) * wing.rootChord;

for i = 1:discretize(2)
    % define planar distribution 
    wing.controlPoint(:, i, 1) = linspace(0, wing.chordDistribution(i), discretize(1));
    wing.controlPoint(:, i, 2) = wing.airfoilCoefficients(1) * sin(1 * pi * linspace(0, wing.chordDistribution(i), discretize(1))) ...
        + wing.airfoilCoefficients(2) * sin(2 * pi * linspace(0, wing.chordDistribution(i), discretize(1)))  ...
        + wing.airfoilCoefficients(3) * sin(3 * pi * linspace(0, wing.chordDistribution(i), discretize(1)));

    wing.controlPoint(:, i, 3) = x(i) * wing.span * cos(wing.dihedral) * cos(wing.sweep);


end

% add third dimension to x and y
wing.controlPoint(:, :, 1) = wing.controlPoint(:, :, 1) + abs(wing.controlPoint(:, :, 3)) * sin(wing.sweep);
wing.controlPoint(:, :, 2) = wing.controlPoint(:, :, 2) + abs(wing.controlPoint(:, :, 3)) * sin(wing.dihedral);

% fancy plots to check wing shape
% scatterControlPoints(wing)
% plotWingSurface(wing)

%% Normal vector 


wing.normal = zeros(discretize(1), discretize(2), 3);
wing.normal(:,:, 1) = ones(discretize(1), discretize(2));

% We may now take advantage of the sine series description of the mean
% line. x-derivative will be obtained analytically.
yPrime = zeros(discretize(2), 1);
for i = 1:discretize(2)
    yPrime(i) = wing.airfoilCoefficients(1) * sin(1 * pi * linspace(0, wing.chordDistribution(i), discretize(1))) ...
        + wing.airfoilCoefficients(2) * sin(2 * pi * linspace(0, wing.chordDistribution(i), discretize(1))) * 2 * pi  ...
        + wing.airfoilCoefficients(3) * sin(3 * pi * linspace(0, wing.chordDistribution(i), discretize(1))) * 3 * pi;
end








%% functions
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


end