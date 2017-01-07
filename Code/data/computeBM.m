function [gHat, sigmaBM] = computeBM(data, fun)

N = size(data, 1);
b = floor(sqrt(N));
a = b;

% compute the general mean
gHat = 0;
for i = 1 : N
    gHat = gHat + fun(data(i, :));
end
gHat = gHat / N;

% compute the means for batches
g = zeros(a, 1);
for k = 1 : a
    for i = 1 : b
        g(k) = g(k) + fun(data((k-1)*b + i, :));
    end
    g(k) = g(k) / b;
end

% Compute the batch means estimator
sigmaBM = 0;
for k = 1 : a
    sigmaBM = sigmaBM + (g(k) - gHat)^2;
end
sigmaBM = sigmaBM * b / (a - 1);