function kappa = buildField_Continuous(theta, D, V, xKL)

nEig = length(theta);
logKappa = 0;
for i = 1 : nEig
   logKappa = logKappa + sqrt(D(i, i)) * theta(i) * V(:, i); 
end

kappa = @(x) interp1(xKL, exp(logKappa), x);


return
