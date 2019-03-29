function F = assembleRHS_quadrature(f, x)

N = length(x) - 1;
F = zeros(N-1, 1);

for j = 2 : N
    FF = @(xx) f(xx) .* (xx-x(j-1))/(x(j) - x(j-1));    
%     xFine = linspace(x(j-1), x(j), 10);
%     F(j-1) = trapz(xFine, FF(xFine));
    F(j-1) = simpson(FF, x(j-1), x(j));  
    FF = @(xx) f(xx) .* (x(j+1)-xx)/(x(j+1) - x(j));
%     xFine = linspace(x(j), x(j+1), 10);
%     F(j-1) = F(j-1) + trapz(xFine, FF(xFine));
    F(j-1) = F(j-1) + simpson(FF, x(j), x(j+1));
end