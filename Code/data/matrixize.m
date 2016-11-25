function A = matrixize(a)

n = length(a);
A = zeros(sqrt(n), sqrt(n));

count = 1;
for i = 1 : sqrt(n)
    for j = 1 : sqrt(n)
        A(j, i) = a(count);
        count = count + 1;
    end
end
