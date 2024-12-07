function gamma = buildLinearSystem(uInf, wing, tail)
% Build linear system to solve wing-tail coupled Weissinger algorithm.

%% 

rCP = reshape(wing.controlPoint, [], 3); 
normal = reshape(wing.normal, [], 3); 
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
VBL = VBL + repmat(uInf', N, 1) * 1e3;
VBR = VBR + repmat(uInf', N, 1) * 1e3;


%% Gamma Matrix
A = zeros(N);
b = zeros(N, 1);
for j = 1:N
    for k = 1:N
        qLeft = q(VFL(j,:) - VBL(j,:), rCP(k,:) - VBL(j,:), rCP(k,:) - VFL(j,:));
        qForward = q(VFL(j,:) - VFR(j,:), rCP(k,:) - VFL(j,:), rCP(k,:) - VFR(j,:));
        qRight = q(VBR(j,:) - VFR(j,:), rCP(k,:) - VFR(j,:), rCP(k,:) - VBR(j,:));

        A(j, k) = dot(qLeft + qForward + qRight, normal(k,:));
    end
end

for k = 1:N
    b(k) = - dot(uInf, normal(k,:));
end

gamma = A\b;






end