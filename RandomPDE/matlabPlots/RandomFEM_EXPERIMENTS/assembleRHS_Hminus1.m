function F = assembleRHS_Hminus1(f1, f2, x)

F = zeros(length(x)-2, 1);
N = length(x) - 1;
for j = 2 : N
    FF = @(xx) f1(xx) .* (xx-x(j-1))/(x(j) - x(j-1));
    F(j-1) = integral(FF, x(j-1), x(j));
    FF = @(xx) f2(xx) * 1/(x(j) - x(j-1));
    F(j-1) = F(j-1) + integral(FF, x(j-1), x(j));
    FF = @(xx) f1(xx) .* (x(j+1)-xx)/(x(j+1) - x(j));
    F(j-1) = F(j-1) + integral(FF, x(j), x(j+1));
    FF = @(xx) f2(xx) * (-1/(x(j+1) - x(j)));
    F(j-1) = F(j-1) + integral(FF, x(j), x(j+1));
end