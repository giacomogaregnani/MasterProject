function Y = EM_1S(Y, f, g, h)

M = length(Y);
xi = randn(M, 1); 
Y = Y + h * f(Y) + sqrt(h) * g(Y) .* xi;

return
