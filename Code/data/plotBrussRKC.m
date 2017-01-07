clc; clear; close all;

hText = [20000, 10000, 5000, 2500, 1250];
h = hText * 1e-5;
nExp = length(h);
meanStages = zeros(size(h));
meanCost = zeros(size(h));
accRatio = zeros(size(h));
sigmaBM = zeros(size(h));
maxStagesTh = zeros(size(h));

% testTwo N = 50; testThree N = 250

x = linspace(0.005, 0.035, 200);
leg = cell(nExp, 1);

figure
hold on
for i = 1 : nExp
   results = dlmread(['Bruss/testTwo_04_01_2017_04_22_', num2str(hText(i)), '.txt']);
   alpha = results(:, 1);
   nStages = results(:, 2);
   accRatio(i) = length(unique(results(:, 1))) / size(results, 1);
     
   meanStages(i) = mean(results(:, 2));
   meanCost(i) = 10 / h(i) * meanStages(i); 
   maxStagesTh(i) = ceil(sqrt(2 * h(i) * (101)^2));
   maxCostTh(i) = 10 / h(i) * maxStagesTh(i);
   
   [f, xi] = ksdensity(alpha, x);
   plot(xi, f);    
   
   leg{i} = ['h = ', num2str(h(i))];
   
   fun = @(x) x;
   [gHat(i), sigmaBM(i)] = computeBM(alpha, fun);
   
   varEstimate(i) = sigmaBM(i) / 50000;   
end

legend(leg)
set(gca, 'yticklabel', '')
xlabel('\alpha')
% 
% figure
% loglog(h, meanCost, 'o-')
% xlim([0.01, 0.3])
