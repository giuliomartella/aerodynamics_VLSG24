clear
clc
close all

% Add the directory containing airfoil data to the path


% Import airfoil profile data
filename = 'naca2412.dat';
airfoil = importXfoilProfile(filename);


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
meanLine.xMl = linspace(0, 1, n); % Uniform x-coordinates from leading to trailing edge.

% Fit cubic splines to the upper and lower surfaces.
ppL = spline(xl, yl); 
ppU = spline(xu, yu); 

% Compute the y-coordinates of the camber line as the average of the surface heights.
ylset = ppval(ppL, meanLine.xMl); 
yuset = ppval(ppU, meanLine.xMl); 
meanLine.yMl = (yuset + ylset) * 0.5; 

nacaC = polyfit(meanLine.xMl, meanLine.yMl, 15);

xPoly = linspace(0, 1);
yPoly = polyval(nacaC, xPoly);
yPrime =  polyval(polyder(nacaC), xPoly);

% Store the original airfoil data in the output structure for reference.
meanLine.x = airfoil.x; 
meanLine.y = airfoil.y; 


% Plot the airfoil, mean camber line, and its derivative
figure
plot(meanLine.x, meanLine.y) % Original airfoil
hold on
axis equal
grid on
plot(meanLine.xMl, meanLine.yMl, 'r') % Mean camber line
plot(xPoly, yPoly, 'g') % Camber line derivative
plot(xPoly, yPrime, 'go') % Camber line derivative


