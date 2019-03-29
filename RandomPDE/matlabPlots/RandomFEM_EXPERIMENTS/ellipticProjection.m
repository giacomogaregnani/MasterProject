function U = ellipticProjection(x, X, u)

A = assembleMatrixProb(X);
AHat = assembleMatrixProj(x, X);
U = [0; A \ (AHat * u(2:end-1)'); 0]';

end