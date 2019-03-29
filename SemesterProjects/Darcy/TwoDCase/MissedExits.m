% Script for explanation of MISSED EXITS 
close all

sigma = 0.8;

V = @(x,y) 2 * [0; x + 1];
dV = @(x,y) 2 * [0; x + 1];
f = @(x,y) -dV(x,y);
g = @(x,y) sigma * eye(2);

X0 = [0.8; 0.8];
Time = [0, 1];
N = [50, 3000];
M = 1;

% plot the SQUARE
p(1) = plot([1,1], [-1,1], 'k', 'LineWidth', 2);
hold on

axis([-1.2 1.2 -1.2 1.2])
axis equal
hold on

W = BrownianMotion2D(Time, N(end), M);

for i = 1:length(N)
    NaiveForPlots(X0, f, g, W(:,1:N(end)/N(i):end), Time, p)
end

legend('Domain', 'Big Timestep', 'Small Timestep', 'Location', 'SouthWest')
% axis([-0.8 1.1 0.5 1])
