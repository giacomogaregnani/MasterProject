D = 10;
accuracy = 10 .^ [-1 : -1 : -8];
rho = 1;

costMC = accuracy.^(-3);
costMLMC = accuracy.^(-2) .* (log2(accuracy)).^2 * D;
costSMLMC = costMLMC .* (1 + sqrt(T * rho) ./ abs(log2(accuracy)));

figure
loglog(accuracy, costMC, '-o')
hold on
loglog(accuracy, costMLMC, '-or')
loglog(accuracy, costSMLMC, '-ok')