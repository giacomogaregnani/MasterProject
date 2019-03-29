function A = assembleMatrixOld(k, x)

N = length(x) - 1;

uDiag = -1 ./ (x(3:end-1) - x(2:end-2));
lDiag = uDiag;
mDiag = 1 ./ (x(3:end) - x(2:end-1)) + 1 ./ (x(2:end-1) - x(1:end-2));

uDiag = [0, uDiag];
lDiag = [lDiag, 0];

A = spdiags([lDiag', mDiag', uDiag'], [-1, 0, 1], N-1, N-1);

end