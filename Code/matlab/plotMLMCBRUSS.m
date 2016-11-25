%% Load numerical solution Monte Carlo h variable N = 100
clear; clc; close all

plotDis = true;

for i = 0 : 2
    resultsfileML{i+1} = ['../data/testBRUSS_sMLMC_15_11_2016_12_21_', ...
        num2str(i), '.txt'];
    costFileML{i+1} = ['../data/testBRUSS_sMLMC_15_11_2016_12_21_', ...
        num2str(i), '_cost', '.txt'];
    resultsfileMC100{i+1} = ['../data/testBRUSS_sMC15_11_2016_12_21_', ...
        num2str(i), '.txt'];
    costFileMC100{i+1} = ['../data/testBRUSS_sMC15_11_2016_12_21_', ...
        num2str(i), '_cost', '.txt'];
    resultsfileMC1000{i+1} = ['../data/testBRUSS_sMC_1000_15_11_2016_12_52_', ...
        num2str(0), '.txt'];
    costFileMC1000{i+1} = ['../data/testBRUSS_sMC_1000_15_11_2016_12_52_', ...
        num2str(0), '_cost', '.txt'];    
    resultsfileMC10{i+1} = ['../data/testBRUSS_sMC_10_15_11_2016_02_08_', ...
        num2str(i), '.txt'];
    costFileMC10{i+1} = ['../data/testBRUSS_sMC_10_15_11_2016_02_08_', ...
        num2str(i), '_cost', '.txt'];   
end
nExperience = length(resultsfileML);

%% Plot distributions

meanMC = zeros(nExperience, 1);
varMC = zeros(nExperience, 1);
meanML = meanMC;
varML = varMC;
meanMC1000 = meanMC;
varMC1000 = varMC;
meanMC10 = meanMC;
varMC10 = varMC;

if plotDis
    figure
    hold on
    for k = 1 : nExperience
        results = dlmread(resultsfileML{k});
        x = results(5001:end);
        meanML(k) = mean(x);
        varML(k) = var(x);
        [f, xi] = ksdensity(x);
        plot(xi, f)
    end
    set(gca, 'YTickLabel', '')
    
    figure
    hold on
    for k = 1 : nExperience
        results = dlmread(resultsfileMC100{k});
        x = results(5001:end);
        meanMC(k) = mean(x);
        varMC(k) = var(x);
        [f, xi] = ksdensity(x);
        plot(xi, f)
    end    
    set(gca, 'YTickLabel', '')  
    
    figure
    hold on
    for k = 1 : nExperience
        results = dlmread(resultsfileMC1000{k});
        x = results(5001:end);
        meanMC1000(k) = mean(x);
        varMC1000(k) = var(x);
        [f, xi] = ksdensity(x);
        plot(xi, f)
    end    
    set(gca, 'YTickLabel', '')
    
    figure
    hold on
    for k = 1 : nExperience
        results = dlmread(resultsfileMC10{k});
        x = results(5001:end);
        meanMC10(k) = mean(x);
        varMC10(k) = var(x);
        [f, xi] = ksdensity(x);
        plot(xi, f)
    end    
    set(gca, 'YTickLabel', '')
end

%% Error - work

for k = 1 : nExperience
    costMC(k) = dlmread(costFileMC100{k});
    costMC(k) = costMC(k) / 50000;
    errMC(k) = abs(meanMC(k) - 1);
    costMC1000(k) = dlmread(costFileMC1000{k});
    costMC1000(k) = costMC1000(k) / 50000;
    errMC1000(k) = abs(meanMC1000(k) - 1);
    costMC10(k) = dlmread(costFileMC10{k});
    costMC10(k) = costMC10(k) / 50000;
    errMC10(k) = abs(meanMC10(k) - 1);
    costML(k) = dlmread(costFileML{k});
    costML(k) = costML(k) / 50000;
    errML(k) = abs(meanML(k) - 1);
end

figure
loglog(errMC, costMC, 'o-')
hold on
loglog(errMC1000, costMC1000, 'o-')
loglog(errMC10, costMC10, 'o-')
loglog(errML, costML, 'o-')
legend('MC', 'MC1000', 'MC10', 'ML')







