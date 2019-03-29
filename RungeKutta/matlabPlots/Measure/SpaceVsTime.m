
results = 0;
nReps = 100;

for i = 1 : nReps
    results = results + dlmread(['Measure1_', num2str(i), '.txt']);
end

results = results / nReps;

ref = 0; %results(1, 2);

diff = abs(results(:, 2) - ref);

figure
loglog(results(:, 1), diff)
hold on
loglog(results(:, 1), 10 * results(:, 1).^(-0.5))
loglog(results(:, 1), 10 * results(:, 1).^(-1))

conv = log2(diff(1:end-1) ./ diff(2:end));

%%
distr = [];
for i = 1 : nReps
   
   distr = [distr; dlmread(['Measure1_', num2str(i), '_full.txt'])];
   
end

figure
histogram(distr, 'normalization', 'pdf')