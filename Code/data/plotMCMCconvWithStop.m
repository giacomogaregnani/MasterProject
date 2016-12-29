clear; clc; close all;

% h = [4000, 2857, 2040, 1470, 1052, 751, 537, 384, 274, 50];
h = [5555, 4000, 2857, 2040, 1470, 1052];
h = h / 1e5;

exactTheta = [0.2, 0.2, 3.0];
gExact = exactTheta * exactTheta';

nExp = length(h);
meanTheta = zeros(nExp, 3);

fileName = 'MeanExp/testHell10_26_12_2016_06_16_'; nExpReps = 1;% MP Fitznag
% fileName = 'MeanExp/testFitznagMP_20_12_2016_09_46_'; nExpReps = 1;% MP Fitznag
% fileName = 'MeanExp/testFitznagEE_14_12_2016_04_35_'; nExpReps = 10;% EE Fitznag


for i = nExp : -1 : 1
    
    display(i)
    
    for k = 1 : nExpReps
        results = dlmread([fileName, num2str(round(h(i)*1e5)), '_', num2str(k-1), '.txt']);
        N = size(results, 1);
        meanTheta(i, :) = mean(results);
        g = sqrt(sum(results.^2, 2));
        
        if i == nExp
            expGThetaExact(k) = mean(g);
            varGThetaExact(k) = var(g);
        else
            meanG(k) = mean(g);
            tmpErrMean(k) = abs(mean(g) - gExact);
            tmpErrVar(k) = abs(var(g) - varExact);
        end
    end
    
    if i == nExp
        gExact = mean(expGThetaExact);
        varExact = mean(varGThetaExact);
    else
        errorGTheta(i) = mean(tmpErrMean);
        errorVar(i) = var(tmpErrMean);
    end
    
    % batch means
    bS = floor(N^(1/3));
    count = 1;
    
    for b = bS
        clear Y
        a = N / b;
        for k = 1 : a
            Y(k) = mean(g((k-1)*b+1 : k*b));
        end
        sigmaHat(i) = b / (a - 1) * sum((Y - mean(g)).^2);
        sigmaHat(i) = sigmaHat(i) / N;
        count = count + 1;
    end    
    
end

figure
loglog(h(1:end-1), errorGTheta, 'o-');
hold on
loglog(h(1:end-1), h(1:end-1).^2, 'k--')
title('Expectation of g(\theta)')
ordMeanG = computeOrder(h(1:end-1), errorGTheta);

% all components
figure
loglog(h(1:end-1), abs(meanTheta(1:end-1, :) - repmat(meanTheta(end, :), nExp-1, 1)), 'o-');
hold on
loglog(h(1:end-1), h(1:end-1), 'k--')

%% Hellinger distance
H = zeros(nExp - 1, 1);

figure
hold on
for k = 1 : 1
    refResults = dlmread([fileName, num2str(round(h(end)*1e5)), '_', num2str(k-1), '.txt']);
    meanRef = mean(refResults);
    SigmaRef = cov(refResults);
    [f, xi] = ksdensity(refResults(:, 1));
    plot(xi, f);
    leg{1} = 'ref';
end

count = 1;
for i = 1 : nExp - 1
    
    display(i)
    
    for k = 1 : nExpReps
        results = dlmread([fileName, num2str(round(h(i) * 1e5)), '_', num2str(k-1), '.txt']);
        Sigma = cov(results);
        Stest{i} = Sigma;
        meanTest = mean(results);
        hTemp(k) = computeNormalHellDistance(meanTest, meanRef, Sigma, SigmaRef);
        if k == 1 && mod(i, 2) == 0
            [f, xi] = ksdensity(results(:, 1));
            plot(xi, f);
            leg{count+1} = ['h = ', num2str(h(i))];
            count = count + 1;
        end
        %         hTemp(k) = computeHellingerDistance(results, refResults, 1e3, false);
    end
    H(i) = mean(hTemp);
    varH(i) = var(hTemp);
    
end
legend(leg)

for i = 1 : nExp - 1
    Coeff(i) = (det(Stest{i}) * det(SigmaRef))^0.25 / sqrt(det(0.5 * (Stest{i} + SigmaRef)));
end
figure
loglog(h(1:end - 1), 1 - Coeff)
hold on
loglog(h(1:end - 1), h(1:end - 1))
% 
% for i = 1 : nExp - 1
% %     ErrVarSingleParam(i, :) = abs(diag(Stest{i}) - diag(SigmaRef));
%     ErrVarSingleParam(i) = norm(Stest{i} - SigmaRef);
% end
% for i = 1 : 3
%     ordSingleVar(i, :) = computeOrder(h(1:end-1)', ErrVarSingleParam(:, i));
% end
% figure
% loglog(h(1:end - 1), ErrVarSingleParam)
% hold on
% loglog(h(1:end-1), h(1:end-1), 'k--')
% loglog(h(1:end-1), h(1:end-1).^2, 'k-')

figure
loglog(h(1:end-1), H, 'o-')
hold on
loglog(h(1:end-1), 10 * h(1:end-1).^1, 'k--')
title('Hellinger')

% Order of convergence
hOrder = computeOrder(h(1:end-1), H');


