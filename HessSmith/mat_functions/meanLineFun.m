function meanLine = meanLineFun(airfoil)
% MEANLINEFUN computes the mean camber line of an airfoil using spline coefficients.

lower = false;
xl = []; yl = []; xu = []; yu = []; 
for i = 1:numel(airfoil.x)
    lower = lower || (airfoil.y(i) < 1e-10);
    if lower
        xl = [xl; airfoil.x(i)];
        yl = [yl; airfoil.y(i)];
    else
        xu = [xu; airfoil.x(i)];
        yu = [yu; airfoil.y(i)];
    end
end

n = 1e2;
meanLine.xMl = linspace(0, 1, n);
ylset = spline(xl, yl, meanLine.xMl);
yuset = spline(xu, yu, meanLine.xMl);
yMl = (yuset + ylset) * 0.5;

% make mean line regular near leading edge
cfit=fit(meanLine.xMl(:),yMl(:),'smoothingspline');
meanLine.yMl = feval(cfit, meanLine.xMl)';

dy = diff(meanLine.yMl) ./ diff(meanLine.xMl);  
dy = [dy(1) dy];  
cfit=fit(meanLine.xMl(:),dy(:),'smoothingspline');
meanLine.dy = feval(cfit, meanLine.xMl)';

% meanLine.dy = dy;

 end
