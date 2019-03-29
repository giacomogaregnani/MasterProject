function phi = ComputeExitProbFD(X0, Time, Bounds, BoundCond, f, sigma)
% Computes the PROBABILITY OF EXIT using finite differences

% Set up of time span and space interval
dx = 0.01;
h = (Time(2)-Time(1))/(2^12);
x = Bounds(1):dx:Bounds(2);
t = Time(1):h:Time(2);

% Define grid
N = length(x) - 2;
M = length(t);

if BoundCond(2) == 0
    
    A = spdiags([[sigma^2/dx * ones(N-1,1) - f(x(3:end-1))';0], ...
        -2 * sigma^2/dx * ones(N,1), ...
        [0;sigma^2/dx * ones(N-1,1) + f(x(2:end-2))']], -1:1, N, N);
    
    A = 1/(2*dx) * A;
    F = 1/(2*dx) * [sigma^2/dx - f(x(2)); zeros(N-2,1); sigma^2/dx + f(x(end-1))];
    u = zeros(N,1);
    theta = 1;
    
    for i = 1 : M - 1
        u = (speye(N) - h * theta * A) \ ((speye(N) + h * (1-theta) * A) * u + h * F);
    end
    u = [1; u; 1];
    
elseif BoundCond(2) == 1
    
    theta = 0.5;
    A = spdiags([[sigma^2/dx * ones(N-1,1) - f(x(3:end-1))'; -1 / (h*theta) ; -1], ...
        [-2 * sigma^2/dx * ones(N,1); 1 / (h * theta)] ...
        [0;sigma^2/dx * ones(N,1) + f(x(2:end-1))']], -1:1, N + 1, N + 1);
    
    A = 1/(2*dx) * A;
    F = 1/(2*dx) * [sigma^2/dx - f(x(2)); zeros(N,1)];
    u = zeros(N + 1,1);
    
    for i = 1 : M - 1
        u = ([speye(N), zeros(N,1); zeros(1,N+1)] - h * theta * A) \ (([speye(N), zeros(N,1); zeros(1,N+1)] + h * (1-theta) * A) * u + h * F);
    end
    
    u = [1 ; u];
end

phi = interp1(x,u,X0);

