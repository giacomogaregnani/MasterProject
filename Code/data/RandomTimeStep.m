clc; clear; close all

h = 0.1 ./ 2 .^ [0 : 7];
hText = floor(h * 1e5);
nExp = length(h);

phi = @(Y) Y;
count = 1;

trueVal = exp(0.5 * 10);

for i = 1 : nExp
    
%     trueVal = exp(0.5 * h(i));
    
    YRK = dlmread(['invMeasRK_hTest', num2str(hText(i)), '.txt']);
    YMP = dlmread(['invMeasMP_hTest', num2str(hText(i)), '.txt']); 
    YEE = dlmread(['invMeasEE_hTest', num2str(hText(i)), '.txt']);

    phiRK = abs(phi(YRK) - trueVal);
    phiMP = abs(phi(YMP) - trueVal);
    phiEE = abs(phi(YEE) - trueVal);
    mPhiRK(count) = mean(phiRK);
    mPhiMP(count) = mean(phiMP);
    mPhiEE(count) = mean(phiEE);
    count = count + 1;
    
end

figure
loglog(h, abs(mPhiRK), 'o-')
hold on
loglog(h, abs(mPhiMP), 'o-')
loglog(h, abs(mPhiEE), 'o-')
% loglog(h, h.^0.5)
loglog(h, h.^0.5 * 80)
loglog(h, h * 100)
loglog(h, h.^1.5 * 80)
loglog(h, h.^2 * 50)
loglog(h, h.^2.5 * 80)
% loglog(h, h.^2.5 * 40)
% loglog(h, h.^3)
legend('RungeKutta', 'MidPoint', 'Euler', '0.5', '1', '1.5', '2', '2.5', 'Location', 'Best')