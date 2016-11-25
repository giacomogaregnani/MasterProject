clear; clc; close all

resultsfile = cell(3, 8);
h = [0.1, 0.05, 0.025];
M = [25, 50, 75, 100, 125, 150, 175, 200];
hText = {'0_1', '0_05', '0_025'};
[nH, nM] = size(resultsfile);

for i = 1 : nH
    for j = 1 : nM
        resultsfile{i,j} = ['../data/testFITZNAG_MC_20_11_2016_03_56__MC_', num2str(M(j)), '_h_', hText{i}, '.txt'];
    end
end

%% Plot distributions

meanMC = cell(nH, nM);
varMC = cell(nH, nM);
errors = cell(nH, nM);

for i = 1 : nH
    figure
    hold on
    leg = cell(nM, 1);
    for j = 1 : nM
        results = dlmread(resultsfile{i, j});
        x = results(5001:end, :);
        meanMC{i, j} = mean(x);
        varMC{i, j} = mean(var(x));
        [f, xi] = ksdensity(x(:, 3));
        plot(xi, f)
        leg{j} = ['M = ', num2str(M(j))];
        errors{i, j} = norm([0.2, 0.2, 3.0] - meanMC{i, j});
    end
    legend(leg);
end
set(gca, 'YTickLabel', '')

% plot error as a fct. of M fixed h
figure
tmpErr = zeros(1, nM);
tmpVar = zeros(1, nM);
orderM = zeros(nH, nM - 1);
for i = 1 : nH
    for j = 1 : nM
        tmpErr(j) = errors{i, j};
        tmpVar(j) = varMC{i, j};
    end
    loglog(M, tmpVar, 'o-');
    hold on
    order(i, :) = - log(tmpVar(2:end)./tmpVar(1:end-1)) ./ log(M(2:end)./M(1:end-1));
end
loglog(M, 0.1./sqrt(M), 'k--');
legend('h = 0.1', 'h = 0.05', 'h = 0.025', 'slope M^{-1/2}')
xlabel('M')

% plot error as a fct. of h fixed M
figure
tmpErrH = zeros(1, nH);
tmpVarH = zeros(1, nH);
orderH = zeros(nM, nH - 1);
for i = 1 : nM
    for j = 1 : nH
        tmpErrH(j) = errors{j, i};
        tmpVarH(j) = varMC{j, i};
    end
    loglog(h, tmpErrH, 'o-');
    hold on
    orderH(i, :) = log(tmpErrH(2:end)./tmpErrH(1:end-1)) ./ log(h(2:end)./h(1:end-1));
end
loglog(h, h, 'k--');









