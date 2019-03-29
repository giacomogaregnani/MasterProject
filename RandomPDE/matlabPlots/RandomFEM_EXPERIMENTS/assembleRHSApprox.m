function F = assembleRHSApprox(f, x, X, p)

N = length(x) - 1;
h = x(2) - x(1);

% Find coefficients
alpha = zeros(1, N+1);
CK = zeros(1, N);
alpha(1) = 0;
for i = 2 : N+1
    alpha(i) = (X(i) - x(i)) / h^p;
    CK(i-1) = alpha(i) - alpha(i-1);
end

F = zeros(length(x)-2, 1);
for j = 2 : N
    
    % Build local map
    B = 1 + CK(j-1) * h^(p-1);
    b = alpha(j-1) * h^p - CK(j-1) * h^(p-1) * x(j-1);
    T = @(x) B*x + b;
    % Integrate
    FF = @(xx) B * f(T(xx)) .* (xx-x(j-1))/(x(j) - x(j-1));
    F(j-1) = integral(FF, x(j-1), x(j));
    
    % Build local map
    B = 1 + CK(j) * h^(p-1);
    b = alpha(j) * h^p - CK(j) * h^(p-1) * x(j);
    T = @(x) B*x + b;
    % Integrate
    FF = @(xx) B * f(T(xx)) .* (x(j+1)-xx)/(x(j+1) - x(j));
    F(j-1) = F(j-1) + integral(FF, x(j), x(j+1));
    
end

