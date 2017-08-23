function K = buildKmatrix(x, y, pF)

K = zeros(length(x), length(y));

for i = 1 : length(x)
    for j = 1 : length(y)
        K(i, j) = pF.k(x(i), y(j));        
    end
end
