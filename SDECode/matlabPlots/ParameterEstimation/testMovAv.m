clc; clear; close all
%% 

h = 0.001;
T = 100;
N = T / h;
tVec = 0:h:T;
M = 1;

W = [zeros(M, 1), cumsum(sqrt(h) * randn(M, N), 2)];
plot(tVec, W');

% Integral mean
Y = cumsum(W, 2) ./ (1:N+1);
hold on
plot(tVec, Y)

% Other formula
YSDE = zeros(M, N+1);
for i = 1 : N
    YSDE(:, i+1) = h/tVec(i+1) * (W(:, i));
end
plot(tVec, YSDE)