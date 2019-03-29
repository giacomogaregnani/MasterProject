clc; clear; close all
%%

yStep = dlmread('testGrowth.txt');
% yProb = dlmread('testGrowth2.txt');
yDet = dlmread('testGrowthdet.txt');

trueH = yDet(1, end);

%%
N = size(yDet, 1);
M = size(yStep, 1) / N;
t = yDet(:, 1);
h = t(2) - t(1);
p = 2.5; q = 2;
d = 4;

% err = zeros(1, N); 
errStep = zeros(1, N);
% err2 = zeros(1, N-1); 
err2Step = zeros(1, N-1);
for i = 1 : M
    idx1 = (i-1) * N + 1;
    idx2 = idx1 + N - 1;
%     trajProb = yProb(idx1:idx2, 2:end);
    trajStep = yStep(idx1:idx2, 2:end);
        
    for j = 2 : N
%        err(j) = err(j) + norm(trajProb(j, 1:end-1) - yDet(j, 2:end-1)); 
       errStep(j) = errStep(j) + norm(trajStep(j, 1:end-1) - yDet(j, 2:end-1));
%        err2(j-1) = err2(j-1) + abs(trajProb(j, end) - trueH);
       err2Step(j-1) = err2Step(j-1) + abs(trajStep(j, end) - trueH);
    end
end
% err = err / M;
% err2 = err2 / M;
errStep = errStep / M;
err2Step = err2Step / M;

%%
% figure
% loglog(t, err)
% hold on
% loglog(t, errStep);
% loglog(t, 10 * h^p * sqrt(t), 'k--')
% loglog(t, 10 * h^p * t.^1.5, 'k')

figure
loglog(t(2:end), abs(yDet(2:end, end) - trueH));
hold on
% loglog(t(2:end), err2);
loglog(t(2:end), err2Step);
hold on
loglog(t, h^(p+q) * sqrt(t))
loglog([t(2), t(end)], [h^q, h^q]) 