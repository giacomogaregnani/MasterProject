function [X, w, likelihood] = PF(IC, f, g, M, T, obs, sigma)

% Initialization
N = length(obs) - 1;
h = T / N;

w = 1/M * ones(1, M);
X = zeros(M, N+1);
X(:, 1) = IC * ones(1, M);

likelihood = 0;

for k = 2 : N+1
    X(:, k) = randsample(X(:, k-1), M, true, w);
    X(:, k) = EM_1S(X(:, k), f, g, h);
    for j = 1 : M
        w(j) = likFunction(X(j, k), obs(k), sigma);
    end 
    likelihood = likelihood + log(mean(w));
    w = w / sum(w);
end


