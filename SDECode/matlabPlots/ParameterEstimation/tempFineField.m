clc; clear; close all
%%

% delta = 100;
sigma = 1;
rho = 0.1;
nu = 1.1;
k = @(t, s) sigma^2 * 2^(1-nu)/gamma(nu) * (sqrt(2*nu) * abs(t-s)/rho)^nu * besselk(nu, sqrt(2*nu) * abs(t-s)/rho);

n = 1000;
xVec = linspace(0, 1, n);

KK = zeros(n);
for i = 1 : n
    for j = 1 : n
        if j == i
            KK(i, j) = sigma^2;
        else
            KK(i, j) = k(xVec(i), xVec(j));
        end
    end
end

figure
hold on
for i = 1 : 5
    
    % [XX, YY] = meshgrid(xVec, xVec);
    
    % KK = k(XX, YY);
    cholK = chol(KK);
    
    randSample = cholK' * randn(n, 1);
    plot(xVec, randSample)
end