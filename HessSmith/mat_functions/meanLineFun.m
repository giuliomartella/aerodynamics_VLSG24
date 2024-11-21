function meanLine = meanLineFun(airfoil)
% MEANLINEFUN computes the mean camber line of an airfoil using spline coefficients.

% correct trailing edge
meanLine.tau = atan2(airfoil.y(1), 1);
rotation = [cos(meanLine.tau), -sin(meanLine.tau); sin(meanLine.tau), cos(meanLine.tau)];
for i = 1:numel(airfoil.x)
    rotated_X = rotation \ [airfoil.x(i); airfoil.y(i)];
    airfoil.x(i) = rotated_X(1) -0.5;  airfoil.y(i) = rotated_X(2); 
end
%airfoil.x = airfoil.x -0.5;

lower = false;
xl = []; yl = []; xu = []; yu = []; 
for i = 1:numel(airfoil.x)
    lower = lower || (airfoil.y(i) < 0);
    if lower
        xl = [xl; airfoil.x(i)];
        yl = [yl; airfoil.y(i)];
    else
        xu = [xu; airfoil.x(i)];
        yu = [yu; airfoil.y(i)];
    end
end


n = 1e4;
meanLine.xMl = linspace(-0.5, 0.5, n);
ppL = spline(xl, yl);
ppU = spline(xu, yu);
ppLPrime = fnder(ppL);
ppUPrime = fnder(ppU);


ylsetPrime = ppval(ppLPrime, meanLine.xMl);
yusetPrime = ppval(ppUPrime, meanLine.xMl);
dy = (yusetPrime + ylsetPrime) * 0.5;


cfit=fit(meanLine.xMl(:),dy(:),'smoothingspline');
meanLine.xMl = linspace(-0.5, 0.5, n*1e2);
meanLine.dy = feval(cfit, meanLine.xMl)';


ylset = ppval(ppL, meanLine.xMl);
yuset = ppval(ppU, meanLine.xMl);
meanLine.yMl = (yuset + ylset) * 0.5;

meanLine.x = airfoil.x;
meanLine.y = airfoil.y;


end
