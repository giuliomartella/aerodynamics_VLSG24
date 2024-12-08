function [wing, tail] = buildLinearSystem(uInf, wing, tail)
% The function buildLinearSystem constructs a linear system to solve the
% coupled wing-tail Weissinger algorithm. It processes the control points,
% vortex locations, and normal vectors for both the wing and the tail 
% (if present), and uses them to build the matrix A and the right-hand side 
% vector b for solving the system. The function computes the aerodynamic
% coefficients gamma and the panel displacements deltaZ for the wing and tail,
% returning them as 2D matrices for each control point. The function supports
% cases where only the wing or both the wing and the tail are involved in the calculation.

%% 

% Reshape control points for the wing into a 2D matrix (one row per control point)
rCP = reshape(wing.controlPoint, [], 3);  

% Reshape the normal vectors at each control point into a 2D matrix
normal = reshape(wing.normal, [], 3);  

% Reshape vortex locations
VFL = reshape(wing.VFL, [], 3); 
VFR = reshape(wing.VFR, [], 3); 
VBL = reshape(wing.VBL, [], 3);  
VBR = reshape(wing.VBR, [], 3); 



if nargin == 2
    tailFlag = false;
elseif nargin == 3
    tailFlag = true;
    tail_rCP = reshape(tail.controlPoint, [], 3);
    tail_normal = reshape(tail.normal, [], 3);
    tail_VFL = reshape(tail.VFL, [], 3);
    tail_VFR = reshape(tail.VFR, [], 3);
    tail_VBL = reshape(tail.VBL, [], 3);
    tail_VBR = reshape(tail.VBR, [], 3);

    % Combine tail and wing data
    rCP = [rCP; tail_rCP];
    normal = [normal; tail_normal];
    VFL = [VFL; tail_VFL];
    VFR = [VFR; tail_VFR];
    VBL = [VBL; tail_VBL];
    VBR = [VBR; tail_VBR];
end

N = size(rCP, 1);
q = @(r0, r1, r2) ...
    (dot(r0, (r2 / norm(r2) - r1 / norm(r1)))) * ...
    (cross(r1, r2) / (4 * pi * norm(cross(r1, r2))^2));



% Fix Backward Vertices
VBL = VBL + repmat(uInf', N, 1) * 1e6;
VBR = VBR + repmat(uInf', N, 1) * 1e6;


%% Gamma Matrix
A = zeros(N);
b = zeros(N, 1);
deltaZ = zeros(N, 1);
for j = 1:N
    for k = 1:N
        % Each horseshoe vortex contributes with three filaments
        qLeft = q(VFL(j,:) - VBL(j,:), rCP(k,:) - VBL(j,:), rCP(k,:) - VFL(j,:));
        qForward = q(VFL(j,:) - VFR(j,:), rCP(k,:) - VFL(j,:), rCP(k,:) - VFR(j,:));
        qRight = q(VBR(j,:) - VFR(j,:), rCP(k,:) - VFR(j,:), rCP(k,:) - VBR(j,:));

        A(j, k) = dot(qLeft + qForward + qRight, normal(k,:));
    end
    deltaZ(j) = norm(VFL(j,:) - VFR(j,:));
end

for k = 1:N
    b(k) = - dot(uInf, normal(k,:));
end

gamma = A\b;



%% Split gamma and deltaZ into wing and tail parts, keeping 2D matrices
if tailFlag
    % Split results into wing and tail sections
    numWingPoints = prod(wing.discretize);  
    
    % Reshape gamma and deltaZ into 2D matrices for wing and tail
    wing.gamma = reshape(gamma(1:numWingPoints), wing.discretize(1), wing.discretize(2));  
    wing.deltaZ = reshape(deltaZ(1:numWingPoints), wing.discretize(1), wing.discretize(2));  
    
    tail.gamma = reshape(gamma(numWingPoints + 1:end), tail.discretize(1), tail.discretize(2));  
    tail.deltaZ = reshape(deltaZ(numWingPoints + 1:end), tail.discretize(1), tail.discretize(2));  
else
    % If no tail, just reshape for wing
    wing.gamma = reshape(gamma, wing.discretize(1), wing.discretize(2));  
    wing.deltaZ = reshape(deltaZ, wing.discretize(1), wing.discretize(2));   
end


end