function M = assembleMassMatrix(x)

N = length(x) - 1;
M = zeros(N-1, N-1);
h = 1 / N;

M = h/3 * eye(N-1) - h/6 * diag(ones(N-2,1), 1) -h/6 * diag(ones(N-2,1), -1);