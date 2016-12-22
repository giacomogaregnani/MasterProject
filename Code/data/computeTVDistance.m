function TV = computeTVDistance(sampleOne, sampleTwo, nEvalPoints)
TV = 1;

for j = 1 : size(sampleOne, 2)
    minVal = min(min(sampleOne(:, j)), min(sampleTwo(:, j)));
    maxVal = max(max(sampleOne(:, j)), max(sampleTwo(:, j)));
    evalPoints = linspace(minVal, maxVal, nEvalPoints);
    [f1(:, j), ~] = ksdensity(sampleOne(:, j), evalPoints);
    [f2(:, j), ~] = ksdensity(sampleTwo(:, j), evalPoints);
    tmp(j) = trapz(evalPoints, abs(f1(:, j) - f2(:, j)));
end

for i = 1 : size(sampleOne, 2)
   TV = TV * tmp(i); 
end

TV = 0.5 * TV;

