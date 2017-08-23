function LK = assembleLK(X0_A, X0_B, X, l)

n = length(X0_A);
m = length(X0_B);
p = length(X);

squaredExp = @(x, xP) exp(-(x - xP)^2 / (2 * l^2));

AK = zeros(n, p);

for i = 1 : n
    for j = 1 : p
        AK(i, j) = squaredExp(X0_A(i), X(j)) * ...
            (-l^2 + (X0_A(i) - X(j))^2) / (l^4);
    end
end

BK = zeros(m, p);

for i = 1 : m
    for j = 1 : p
        BK(i, j) = squaredExp(X0_B(i), X(j));
    end
end

LK = [AK; BK];