function wing = postprocessing(wing, tail)

if nargin == 1
    wing = take(wing);
else
    wing = take(wing);
    tail = take(tail);
    wing.cL = (wing.cL * wing.S + tail.cL * tail.S) / (wing.S + tail.S);
    wing.cD = (wing.cD * wing.S + tail.cD * tail.S) / (wing.S + tail.S);
end
    function wing = take(wing)
        wing.gammaDistribution = sum(wing.gamma, 1);
        wing.L2D = wing.gammaDistribution * cos(wing.dihedral) * 1.225;
        wing.D2D = wing.L2D .* sin(- wing.alphaI);
        % wing.L2D = wing.L2D .* cos(wing.alphaI);

        wing.L = sum(wing.L2D.* wing.span ./ wing.discretize(2));
        wing.D = sum(wing.D2D.* wing.span ./ wing.discretize(2));

        wing.L = sum(wing.L2D.* wing.deltaZ(1, :));
        wing.D = sum(wing.D2D.* wing.deltaZ(1, :));

        wing.cL = 2 * wing.L / wing.S / 1.225;
        wing.cD = 2 * wing.D / wing.S / 1.225;


    end

% function wing = take(wing)
%       wing.cL2D = zeros(wing.discretize(2), 1);
% 
%       for j = 1:wing.discretize(2)
%           wing.gammaDistribution(j) = sum(wing.gamma(:,j)) / wing.chordDistribution(j);
%       end
% 
%       wing.cL2D = wing.gammaDistribution * cos(wing.dihedral);
%       wing.cD2D = abs(wing.cL2D .* sin(wing.alphaI));
%       wing.cL2D = wing.cL2D .* cos(wing.alphaI);
% 
%       wing.cL = sum(wing.cL2D) / wing.span;
%       wing.cD = sum(wing.cD2D) / wing.span;
%   end

end