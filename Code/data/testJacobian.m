%% Jacobian of the matrix function

data = dlmread('jacobianHires.txt');
J = data(1 : 8, :);
initialVar = data(9 : end, :);
I = eye(8);

varPrime = @(P) kron(I, J) * P(:) + kron(I, P) * J(:) + 0.01 * I(:);

w = varPrime(initialVar);
for i = 1 : 1000
    v = w / norm(w);
    w = varPrime(matrixize(v));   
    minEig = v' * w / (v' * v); 
end

