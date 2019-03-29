clear
h = 0.005;
N = 1 / h;

for j = 1 : N-1
    lambda(j) = sin(pi * j / (2 * N));
    
    for i = 1 : N-1
        V(j, i) = sin(i*j*pi / N);
    end
    
end

lambda = lambda * (-4 / h^2);
V = V * sqrt(2 / N);

V = [zeros(1, size(V, 2));
     V;
     zeros(1, size(V, 2))];

x = 0 : h : 1;
plot(x, V(:, 1:5))


%% KL 
lambdaINV = 1 ./ sqrt(-lambda);

coeffs  = randn(N-1, 1);
func = 0;
for i = 1 : N-1
   func = func + coeffs(i) * lambdaINV(i) * V(:, i); 
end

figure
plot(x, func)

