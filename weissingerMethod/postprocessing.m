function wing = postprocessing(wing)

wing.cL2D = zeros(wing.discretize(2), 1);

for j = 1:wing.discretize(2)
    wing.cL2D(j) = sum(wing.gamma(:,j)) / wing.chordDistribution(j) * cos(wing.dihedral);
end

wing.cL = sum(wing.cL2D) / wing.span;



end