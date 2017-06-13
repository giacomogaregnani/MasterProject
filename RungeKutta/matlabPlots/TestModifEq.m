clc; clear; close all;

%% Time parameters

T = 30;
h = 0.01;
N = T / h;
timeVec = linspace(0, T, N + 1);

%% Function (space) definition

sigma = 10; rho = 28; theta = 8/3;
f = @(y) [sigma * (y(2) - y(1));
    y(1) * (rho - y(3)) - y(2);
    y(1) * y(2) - theta * y(3)];
initCond = [-10; -1; 40];

% a = 0.2; b = 0.2; c = 3.0;
% f = @(y) [c * (y(1) - y(1)^3 / 3 + y(2));
%           -1/c * (y(1) - a + b * y(2))];
% initCond = [-1; 1];

fAuto = @(t, y) f(y);

%% Monte carlo (modified equation)

figure 
hold on

for i = 1 : 1;
    
    % Simulate theta
    Theta = 2 * rand(N, 1);  
    thetaFunc = @(t) Theta(round(t / h) + 1);
    
    % Build modified function
    F = @(t, y) thetaFunc(t) * f(y);
    
    % Initialize solution
    Y = zeros(length(initCond), N+1);
    Y(:, 1) = initCond;
    
    % Integrate in time (RTS)
    for k = 1 : N-1
        Y(:, k+1) = rk4(fAuto, h * Theta(k), timeVec(k), Y(:, k));
    end
   
    plot(timeVec, Y(1, :), 'color', [0.4, 0.4, 0.4])   
    
end
