clc; clear; close all

h = 0.1 ./ 2 .^ [0 : 11];
hText = floor(h * 1e5);

phi = @(Y) Y;
% phi = @(Y) sum(Y.^2, 2);
count = 1;

trueVal = exp(0.5 * 10);

for H = hText
    
   YRK = dlmread(['invMeasRK_hTest', num2str(H), '.txt']);
   YEE = dlmread(['invMeasEE_hTest', num2str(H), '.txt']);
%    figure 
%    scatter3(YEE(:, 1), YEE(:, 2), YEE(:, 3), '.')
%    hold on
%    scatter3(YRK(:, 1), YRK(:, 2), YRK(:, 3), '.')    
   
   phiRK = phi(YRK);
   phiEE = phi(YEE);
   mPhiRK(count) = mean(phiRK);
   mPhiEE(count) = mean(phiEE);
   count = count + 1;
      
end

figure
loglog(h, abs(mPhiRK - trueVal), 'o-')
hold on
loglog(h, abs(mPhiEE - trueVal), 'o-')
loglog(h, h, 'k--')
loglog(h, h.^2, 'k')
legend('MidPoint', 'Euler', 'Location', 'best')