function LbK = assembleLbK(X0_A, X0_B, X, l)

n = length(X0_A);
m = length(X0_B);
p = length(X);

AbK = zeros(p, n);

squaredExp = @(x, xP) exp(-(x - xP)^2 / (2 * l^2));

for i = 1 : p
    for j = 1 : n
        AbK(i, j) = squaredExp(X(i), X0_A(j)) * ...
            (-l^2 + (X(i) - X0_A(j))^2) / (l^4);
    end
end

BbK = zeros(p, m);

for i = 1 : p
    for j = 1 : m
        BbK(i, j) = squaredExp(X(i), X0_B(j));
    end
end

LbK = [AbK, BbK];
