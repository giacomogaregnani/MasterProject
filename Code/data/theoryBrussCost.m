N = 500;
h = 10 .^ [-1, -3, -5, -7]
alpha = 1;
T = 10;

% Explicit Euler
hEEmin = 1 / (8 * alpha * (N + 1)^2);
costEEmin = T / hEEmin;
costEE = max(costEEmin, T./h);

% RungeKuttaCheb
s = max(2, ceil(sqrt(0.5 * h * 4 * alpha * (N + 1)^2)));
costRKC = T * s ./ h;
