function F = assembleRHSTilde(f, X, x)

F = zeros(length(X)-2, 1);
N = length(X) - 1;
for j = 2 : N
    phi = @(xx) interp1([-1, X(j-1), X(j), X(j+1), 2], [0, 0, 1, 0, 0], xx);
    Phi = phi([x(j-1), x(j), x(j+1)]);
    rPhi = @(xx) interp1([-1, x(j-1), x(j), x(j+1), 2], [0, Phi, 0], xx);
    FF = @(xx) f(xx) .* rPhi(xx);
    if j > 1
        F(j-1) = integral(FF, X(j-1), X(j));
    end
    if j < N+1
        F(j-1) = F(j-1) + integral(FF, X(j), X(j+1));
    end
end