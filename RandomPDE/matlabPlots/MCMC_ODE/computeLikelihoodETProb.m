function L = computeLikelihoodETProb(theta, tObs, yObs, sigmaObs, h, p, nMC)

L = ones(nMC, 1);

for k = 1 : nMC
    
    for i = 1 : length(tObs)
        
        N = round(tObs(i) / h);
        
        y = 1;
        for j = 1 : N
            H = h + h^p * (-1 + 2*rand(1));
            y = y * (1 + H * theta + H^2 * theta^2 / 2);
        end
        L(k) = L(k) * exp(-1 / (2 * sigmaObs^2) * (yObs(i) - y)^2);
        
    end
    
end

L = mean(L);

% L = 1;
% 
% for i = 1 : length(tObs)
%     y = 1 + h * theta + h^2 * theta^2 / 2 + 1 / 3 * h^(2*p);
%     L = L * exp(-1 / (2 * sigmaObs^2) * (yObs(i) - y^(tObs / h))^2);    
% end