function V = normalKernelDensity(data, x)
% V = normalKernelDensity(data, x)
% Each row of data has to be a d-dimensional observation
% Each row of x is a point where the density has to be evaluated
% V is the evaluated density

d = size(data, 2);
n = size(data, 1);
nEvalPoints = size(x, 1);

% Estimate bandwidth (Silverman's rule of thumb)
h = (4 / (n * (d + 2))) ^ (1 / (d + 4))^2 * var(data);
detH = prod(h);

K = @(X) (2 * pi)^(-d / 2) * detH^(-d / 2) * exp(-0.5 * sum(X.^2 ./ repmat(h, n, 1), 2));

V = zeros(nEvalPoints, 1);
parfor j = 1 : nEvalPoints
    V(j) = sum(K(repmat(x(j, :), n, 1) - data));
end
V = V / n;

return