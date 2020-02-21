function f = driftKL(x, alpha, n)

f = 0;
for i = 1 : n
    f = f + alpha(i) * x.^(2*i);
end