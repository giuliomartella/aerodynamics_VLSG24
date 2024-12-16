function wing = postprocessing(wing)

wing.cL2D = zeros(wing.discretize(2), 1);

for j = 1:wing.discretize(2)
    wing.gammaDistribution(j) = sum(wing.gamma(:,j)) / wing.chordDistribution(j);
end

wing.cL2D = wing.gammaDistribution * cos(wing.dihedral);
wing.cD2D = wing.cL2D .* sin(- wing.alphaI);
wing.cL2D = wing.cL2D .* cos(wing.alphaI);

wing.cL = sum(wing.cL2D) / wing.span;
wing.cD = sum(wing.cD2D) / wing.span;



end