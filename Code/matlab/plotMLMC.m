%% Comparison between MC and MLMC
clear
close all
results = dlmread('../data/errors1.txt');
for i = 2 : 8  
    results = results + dlmread(['../data/errors' num2str(i) '.txt']);
end
results = results / 8;
[n, ~] = size(results);

resultsMLMC = results(1 : n/2, 2:-1:1);
resultsMC = results(n/2 + 1 : end, 2:-1:1);

loglog(resultsMLMC(:, 1), resultsMLMC(:, 2),'o-')
hold on
loglog(resultsMC(:, 1), resultsMC(:, 2), 'or-')
legend('MLMC', 'MC');
xlabel('Cost')
ylabel('Accuracy')

epsilon = resultsMLMC(:, 2);
theoryMLMC = abs(log2(epsilon.^0.25)) .* epsilon.^(-2);
theoryMC = epsilon.^(-3);
%loglog(theoryMLMC, epsilon);
%loglog(theoryMC, epsilon);
