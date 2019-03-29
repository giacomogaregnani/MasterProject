function p = ComputeExitProbability(X1, X2, sigma, h)
% Compute EXIT PROBABILITY between TWO CONSECUTIVE TIMESTEPS

% count walls 1 left, 2 bottom, 3 right, 4 up

fp = @(n,z,x1,x2) exp(-2 * (n' * (x1 - z)) * (n' * (x2 - z)) / (h * sigma^2));

n = [-1;0];
z = [-1;X1(2)];
p(1) = fp(n,z,X1,X2);

n = [0;-1];
z = [X1(1);-1];
p(2) = fp(n,z,X1,X2);

n = [1;0];
z = [1;X1(2)];
p(3) = fp(n,z,X1,X2);

n = [0;1];
z = [X1(1);1];
p(4) = fp(n,z,X1,X2);


end