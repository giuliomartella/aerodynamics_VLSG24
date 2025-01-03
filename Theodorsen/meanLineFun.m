function meanLine = meanLineFun(airfoil)
% MEANLINEFUN computes the mean camber line of an airfoil using spline interpolation.
% Additionally, it adjusts the airfoil's trailing edge to ensure proper alignment 
% and calculates the derivative of the camber line for further aerodynamic analysis.

% Correct the trailing edge alignment by calculating the trailing edge angle (tau).
meanLine.tau = atan2(airfoil.y(1), 1); % Calculate the trailing edge angle relative to the chord line.
rotation = [cos(meanLine.tau), -sin(meanLine.tau); sin(meanLine.tau), cos(meanLine.tau)]; 
% Construct a rotation matrix to align the trailing edge horizontally.

for i = 1:numel(airfoil.x)
    rotated_X = rotation \ [airfoil.x(i); airfoil.y(i)]; 
    airfoil.x(i) = rotated_X(1) - 0.5; % Translate x to center the airfoil around x = 0.
    airfoil.y(i) = rotated_X(2); 
end

% Split the airfoil data into upper and lower surfaces based on y values.
lower = false; % Logical flag to indicate transition to the lower surface.
xl = []; yl = []; % Arrays to store lower surface points.
xu = []; yu = []; % Arrays to store upper surface points.

for i = 1:numel(airfoil.x)
    lower = lower || (airfoil.y(i) < 0); % Once y < 0, assume we are on the lower surface.
    if lower
        xl = [xl; airfoil.x(i)]; 
        yl = [yl; airfoil.y(i)]; 
    else
        xu = [xu; airfoil.x(i)]; 
        yu = [yu; airfoil.y(i)];
    end
end

% Define a high-resolution x-coordinate grid for the camber line.
n = 5e4; % Number of points for the camber line.
meanLine.xMl = linspace(-0.5, 0.5, n); % Uniform x-coordinates from leading to trailing edge.

% Fit cubic splines to the upper and lower surfaces.
ppL = spline(xl, yl); 
ppU = spline(xu, yu); 

% Calculate the derivatives of the upper and lower surface splines.
ppLPrime = fnder(ppL); 
ppUPrime = fnder(ppU); 

% Compute the derivative of the camber line as the average of the surface derivatives.
% Computing the average of the derivatives has lower error than numerically
% derive the camber line points.
ylsetPrime = ppval(ppLPrime, meanLine.xMl); 
yusetPrime = ppval(ppUPrime, meanLine.xMl); 
dy = (yusetPrime + ylsetPrime) * 0.5; % Average derivative using derivation linearity.

% Smooth the camber line derivative using a smoothing spline fit.
cfit = fit(meanLine.xMl(:), dy(:), 'smoothingspline'); 
meanLine.xMl = linspace(-0.5, 0.5, n*1e2); % Increase resolution for smoother results.
meanLine.dy = feval(cfit, meanLine.xMl)'; 

% Compute the y-coordinates of the camber line as the average of the surface heights.
ylset = ppval(ppL, meanLine.xMl); 
yuset = ppval(ppU, meanLine.xMl); 
meanLine.yMl = (yuset + ylset) * 0.5; 

% Store the original airfoil data in the output structure for reference.
meanLine.x = airfoil.x; 
meanLine.y = airfoil.y; 

end
