function LLbK = assembleLLbK(X0_A, X0_B, l)

n = length(X0_A);
m = length(X0_B);

AAbK = zeros(n, n);

squaredExp = @(x, xP) exp(-(x - xP)^2 / (2 * l^2));

for i = 1 : n
    for j = 1 : n
        AAbK(i, j) = squaredExp(X0_A(i), X0_A(j)) * ...
            (3 / (l^4) - 6 * (X0_A(i) - X0_A(j))^2 / (l^6) + ...
             (X0_A(i) - X0_A(j))^4 / (l^8));
    end
end

ABbK = zeros(n, m);

for i = 1 : n
    for j = 1 : m
        ABbK(i, j) = squaredExp(X0_A(i), X0_B(j)) * ...
            (-l^2 + (X0_A(i) - X0_B(j))^2) / (l^4);
    end
end

BAbK = zeros(m, n);

for i = 1 : m
    for j = 1 : n
        BAbK(i, j) = squaredExp(X0_B(i), X0_A(j)) * ...
            (-l^2 + (X0_B(i) - X0_A(j))^2) / (l^4);
    end
end

BBbK = zeros(m, m);

for i = 1 : m
    for j = 1 : m
        BBbK(i, j) = squaredExp(X0_B(i), X0_B(j));
    end
end

LLbK = [AAbK, ABbK; BAbK, BBbK];


end
