function A = assembleMatrixProj(x, X)

N = length(x) - 1;
A = zeros(N-1, N-1);

xFine = linspace(0, 1, 50*N);

for i = 2 : N
    for j = 2 : N
        gradPhiI = (xFine > x(i-1)) .* (xFine < x(i)) .* 1 / (x(i) - x(i-1)) + ...
            (xFine >= x(i)) .* (xFine < x(i+1)) .* (-1) / (x(i+1) - x(i));
        gradPhiJ = (xFine > X(j-1)) .* (xFine < X(j)) .* 1 / (X(j) - X(j-1)) + ...
            (xFine >= X(j)) .* (xFine < X(j+1)) .* (-1) / (X(j+1) - X(j));
        A(i-1,j-1) = trapz(xFine, gradPhiI .* gradPhiJ);
    end
end

end