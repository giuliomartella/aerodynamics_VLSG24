function meanLine = meanLineFun(airfoil)
% MEANLINEFUN computes the mean camber line of an airfoil using spline coefficients.
% The function takes an airfoil structure with fields x (chordwise positions) 
% and y (corresponding vertical positions) as input and returns the mean camber line.
% YOU MAY NEED TO ADJUST COEFFICIENTS FOR BAD FORMATTED AIRFOILS .DAT FILES
% IF POINTS ARE SPLITTED IN UPPER AND LOWER SEGMENTS WE SUGGEST TO USE
% MEANCAMBERLINE

% Extract airfoil coordinates
xA = airfoil.x; % x-coordinates of the airfoil
yA = airfoil.y; % y-coordinates of the airfoil2

% Define the evaluation points along the chord
x = linspace(0, 1, 100); 
ySpan = linspace(min(yA), max(yA), 50); % 50 equally spaced points along the y-axis span

% Initialize the output array
y = zeros(size(x)); % Mean line y-coordinates corresponding to x

% Loop through each chordwise point to compute the mean line
for i = 1:100
    % Identify points on the airfoil close to the current x position
    % gap = 5e-3;
    % condition = (xA > (x(i) - 5e-2)) & (xA < (x(i) + 5e-2)); % Logical vector

    gap = 1e-3; % Initial gap size
    
    % Dynamically adjust the gap until condition is not all zeros
    condition = (xA > (x(i) - gap)) & (xA < (x(i) + gap)); % Logical vector
    stat = nnz(condition) > 10;
    while stat == false
        % Identify points on the airfoil close to the current x position
        gap = gap * 1.5; % Increase the gap size (e.g., by 50%)
        condition = (xA > (x(i) - gap)) & (xA < (x(i) + gap)); % Logical vector
        stat = nnz(condition) > 10;
    end
    
    % Extract the nearby points
    points = [xA(condition) yA(condition)];
    np = size(points, 1); % Number of nearby points
    
    % Compute the distance of each ySpan point from the airfoil points
    yDistance = zeros(50, np); % Initialize distance matrix
    for j = 1:50
        for k = 1:np
            yDistance(j, k) = abs(ySpan(j) - points(k, 2)); % Vertical distance
        end
    end
    
    % Find the ySpan point that maximizes the minimum distance to the airfoil points
    [~, index] = min(max(yDistance, [], 2)); 
    y(i) = ySpan(index); % Store the corresponding y value
end

xS = linspace(0, 0.99, 1e4);

cfit=fit(x(:),y(:),'smoothingspline');
yS = feval(cfit, xS)';

yDer = diff(yS) ./ diff(xS);  
yDer = [yDer(1) yDer];  

meanLine = [xS; yS; yDer]'; 




end



