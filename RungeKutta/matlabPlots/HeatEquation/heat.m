clc; clear; close all;

%% Heat equation - 1D FD 

dx = 0.01;
T = 1;
x = 0 : dx : 1;
N = length(x) - 1;

A = -gallery('tridiag', N - 1);
% A(1, end) = 1;
% A(end, 1) = 1;
A = A / (dx^2);

% Eigen functions of the Laplacian
% nEig = 40;
% [V, D] = laplEig(x, nEig);
% D = -1 ./ D;
% sqrtD = diag(sqrt(D));

% u0 = sin(2 * pi * x)';
% u0 = exp(-(x-0.5).^2 / (2 * 0.01))';
u0 = (x < 0.9) .* (x > 0.7) + (x < 0.3) .* (x > 0.1);
plot(x, u0)

%%
f = @(t, x) 1e-2 * A * x;

[t, y] = ode23s(f, [0, T], u0(2:end-1));

y = [zeros(size(t)), y, zeros(size(t))];

figure
% for i = 1 : length(t)
%    plot(x, y(i, :))
%    if i == 1
%        ylim = get(gca, 'ylim');
%    end
%    set(gca, 'ylim', ylim);
%    pause(0.01)
% end

plot(x, y(end, :));