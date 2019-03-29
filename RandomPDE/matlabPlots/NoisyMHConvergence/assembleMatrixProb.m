function A = assembleMatrixProb(x)

N = length(x) - 1;

upperDiag = [0, -1 ./ (x(3:end-1) - x(2:end-2))];
lowerDiag = [-1 ./ (x(3:end-1) - x(2:end-2)), 0];
diag = 1 ./ (x(2:end-1) - x(1:end-2)) + 1 ./ (x(3:end) - x(2:end-1));

A = spdiags([lowerDiag', diag', upperDiag'], [-1, 0, 1], N-1, N-1);

end