function [U] = vortexInfluence(ControlPoint, Extreme_1, Extreme_2, print)


if(nargin == 3)
    print = false;
end


toll = 1e-10;


r0 = Extreme_2 - Extreme_1;
r1 = ControlPoint - Extreme_1;
r2 = ControlPoint - Extreme_2;

CP = cross(r1, r2);

r1_norm = r1 / norm(r1);
r2_norm = r2 / norm(r2);


CP_norm = dot(CP, CP);

if(CP_norm < toll)
    CP_norm = toll;
end

dot_result = dot(r0, r1_norm - r2_norm);

U_Squared = (1/(4*pi)) * dot(r0, r1_norm - r2_norm) .* CP ./ (CP_norm);
% U_NotSquared = (1/(4*pi)) * dot(r0, r1_norm - r2_norm) .* CP ./ (sqrt(CP_norm));

if(print)
%   r1  
%   r2
%   r0
%   r1_norm
%   r2_norm
%   CP_norm
%   dot_result
%   
%   CP ./ (sqrt(CP_norm))
%   dot_result .* CP ./ (sqrt(CP_norm))
end

U = U_Squared;