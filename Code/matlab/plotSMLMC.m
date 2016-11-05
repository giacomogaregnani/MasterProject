%% Comparison between MC and MLMC
clear
close all

results = dlmread('../data/resultsSMLMC1.txt');
for i = 1 : 3
    newR = dlmread(['../data/resultsSMLMC' num2str(i+1) '.txt']);
    results = results + newR;
end
results = results / 4;
[n, ~] = size(results);

resultsMLMC = results(1 : n, 2:-1:1);

loglog(resultsMLMC(:, 1), resultsMLMC(:, 2),'o-')
xlabel('Cost')
ylabel('Accuracy')

% firstResultsSMCMC for problem y' = -100 * y for 0 \leq t \leq 10