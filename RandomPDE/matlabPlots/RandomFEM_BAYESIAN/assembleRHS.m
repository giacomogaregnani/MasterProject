function F = assembleRHS(f, x)

F = zeros(length(x)-2, 1);
N = length(x) - 1;
for j = 2 : N
    FF = @(xx) f(xx) .* (xx-x(j-1))/(x(j) - x(j-1));
    F(j-1) = integral(FF, x(j-1), x(j));
    FF = @(xx) f(xx) .* (x(j+1)-xx)/(x(j+1) - x(j));
    F(j-1) = F(j-1) + integral(FF, x(j), x(j+1));
end