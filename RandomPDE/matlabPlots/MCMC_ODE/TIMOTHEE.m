h = 1;
lambda = -3;
mu = sqrt(2);

X0 = 1;

for i = 1 : 10000
    X(i) = 1 / (1 - lambda * h - mu * sqrt(h) * randn) * X0;
end

EX = mean(abs(X));

plot(abs(X))

Xt = repmat(X0, 1, 10000);

for j = 1 : 10000
    for i = 1 : 1
        Xt(i+1, j) = 1 / (1 - lambda * h - mu * sqrt(h) * randn) * Xt(i, j);
    end
end

figure
plot(abs(Xt))