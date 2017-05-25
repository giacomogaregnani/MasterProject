%% Solution in FD of PDE u_t + u u_x + u_xx + u_xxxx = 0

clc; clear; close all

%% Space Discretization Parameters$

L = 22;
spaceSpan = [-L/2, L/2];

h = 0.1;
N = L / h;

%% Build the matrices 

% Second derivative
B = -gallery('tridiag', N);

% Fourth derivative
c1 = [ones(N-2, 1); 0; 0];
c2 = -4 * [ones(N-1, 1); 0];
c3 = 6 * ones(N, 1);
c4 = -4 * [0; ones(N-1, 1)];
c5 = [0; 0; ones(N-2, 1)];
C = spdiags([c1, c2, c3, c4, c5], [-2 : 2], N, N);
clear c1 c2 c3 c4 c5

% Nonlinear term
a1 = -[ones(N-1, 1); 0];
a2 = [0; ones(N-1, 1)];
A = spdiags([a1, a2], [-1, 1], N, N);
clear a1 a2

%% Impose periodic boundary conditions

% Second derivative
B(1, end) = 1;
B(end, 1) = 1;
B = B / h^2;

% Fourth derivative
C(1, end-1) = 1;
C(1, end) = -4;
C(2, end) = 1;
C(end-1, 1) = 1;
C(end, 1) = -4;
C(end, 2) = 1;
C = C / h^4;

% Nonlinear term
A(1, end) = -1;
A(end, 1) = 1;
A = @(x) x .* (A * x) / (2 * h);

%% Initial condition and solve

xValues = [spaceSpan(1) : h : spaceSpan(2)]';
uZero = sin(2 * pi * xValues(1 : end - 1) / L);
% uZero = zeros(length(xValues) - 1, 1); uZero(2 : 5) = 0.6;

timeSpan = [0, 200];
fODE = @(t, x) -1 * (A(x) + B * x + C * x);

[T, U] = ode23s(fODE, timeSpan, uZero);

% Apply periodic boundaries
U = [U, U(:, 1)];

%% Plot

contourf(xValues, T, U, 20)
colorbar

%% Animation

figure
for i = 1 : 50 : length(T)
   plot(xValues, U(i, :))
   axis([-L/2, L/2, -3, 3])
   title(['T = ', num2str(T(i))])
   pause(0.0001)    
end