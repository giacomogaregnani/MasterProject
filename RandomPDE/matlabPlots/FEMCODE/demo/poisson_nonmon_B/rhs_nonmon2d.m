function [ out ] = rhs_nonmon2d(x, der)
%% RHS FUNCTION FOR NONLINEAR NONMONOTONE PROBLEM

if (nargin == 1) || (0 == sum(der))
    out = -0.8e1 .* (sin(0.8e1 .* pi .* sin(pi .* x(:,1)) .* x(:,2) .* (1 - x(:,2)))...
    + 0.8e1 .* x(:,1) .* cos(0.8e1 .* pi .* sin(pi .* x(:,1)) .* x(:,2) .* (1 - x(:,2)))...
    .* pi .^ 2 .* cos(pi .* x(:,1)) .* x(:,2) .* (1 - x(:,2))) .* cos(pi .* x(:,1))...
    .* pi .* x(:,2) .* (1 - x(:,2)) + 0.8e1 .* (0.1e1 + x(:,1) .* sin(0.8e1 .* pi...
    .* sin(pi .* x(:,1)) .* x(:,2) .* (1 - x(:,2)))) .* pi .^ 2 .* sin(pi .* x(:,1))...
    .* x(:,2) .* (1 - x(:,2)) - (0.8e1 .* sin(pi .* x(:,1)) .* (1 - x(:,2))...
    - 0.8e1 .* sin(pi .* x(:,1)) .* x(:,2)) .^ 2 ./ (0.1e1 + 0.64e2...
    .* sin(pi .* x(:,1)) .^ 2 .* (x(:,2) .^ 2) .* ((1 - x(:,2)) .^ 2)) + 0.16e2...
    .* (0.2e1 + atan(0.8e1 .* sin(pi .* x(:,1)) .* x(:,2) .* (1 - x(:,2))))...
    .* sin(pi .* x(:,1));
else 
    out = zeros(size(x,1), 1);
end

end