function H = computeHellingerDistance(sample1, sample2, nEvalPoints, plotResults)

f1 = zeros(nEvalPoints, size(sample1, 2));
f2 = f1;
evalPoints = zeros(nEvalPoints, size(sample1, 2));

for j = 1 : size(sample1, 2)
    minVal = min(min(sample1(:, j)), min(sample2(:, j)));
    maxVal = max(max(sample1(:, j)), max(sample2(:, j)));
    evalPoints(:, j) = linspace(minVal, maxVal, nEvalPoints);
    [f1(:, j), ~] = ksdensity(sample1(:, j), evalPoints(:, j));
    [f2(:, j), ~] = ksdensity(sample2(:, j), evalPoints(:, j));    
    if plotResults && j == 2
        figure
        hold on
        plot(evalPoints, f1(:, j))
        plot(evalPoints, f2(:, j))
        xlim([2, 3.5])
    end
end

% Compute Hellinger distance
integral = 1;
for j = 1 : size(sample1, 2)
    integral = integral * trapz(evalPoints(:, j), (sqrt(f1(:, j)) - sqrt(f2(:, j))).^2);
end
H = sqrt(0.5 * integral);

end

