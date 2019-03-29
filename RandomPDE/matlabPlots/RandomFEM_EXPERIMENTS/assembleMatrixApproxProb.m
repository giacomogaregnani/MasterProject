function A = assembleMatrixApproxProb(x, X, p)

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

% Diagonal
A =  1/h * spdiags((1./(1+h^(p-1)*CK(1:end-1)) + 1./(1+h^(p-1)*CK(2:end)))', 0, N-1, N-1);
% Upper diagonal
A = A + (-1/h) * spdiags([0; 1./(1+h^(p-1)*CK(2:end))'],  1, N-1, N-1);
% Lower diagonal
A = A + (-1/h) * spdiags([1./(1+h^(p-1)*CK(2:end))'; 0], -1, N-1, N-1);

end