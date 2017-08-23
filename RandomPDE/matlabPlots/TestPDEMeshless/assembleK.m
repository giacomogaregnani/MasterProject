function K = assembleK(X1, X2, l)

n = length(X1);
m = length(X2);

K = zeros(n, m);

squaredExp = @(x, xP) exp(-(x - xP)^2 / (2 * l^2));

for i = 1 : n
    for j = 1 : m
        K(i, j) = squaredExp(X1(i), X2(j));
    end
end