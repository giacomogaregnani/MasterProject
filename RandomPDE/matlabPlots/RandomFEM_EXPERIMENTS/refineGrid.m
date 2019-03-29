function [xNew, nFine, nCoarse] = refineGrid(x, err, tol, solNorm)

mark = []; coarse = [];
newPoints = [];
xNew = x; xTmp = x;
N = length(x)-1;

for i = 1 : N
    if err(i) > 1.5 * tol * solNorm / sqrt(N)
        mark = [mark i];
    elseif err(i) < 0.25 * tol * solNorm / sqrt(N)
        coarse = [coarse i];
    end
end
for j = mark
    newPoints = [newPoints (x(j+1) + x(j)) / 2];
end

idx = coarse(diff(coarse)==1)+1;
xTmp(idx) = []; % Take out a point if two consecutive elements have to be coarsened
xNew = sort([xTmp, unique(newPoints)]);

nFine = length(newPoints);
nCoarse = length(idx);

return