% Read results

clc; clear;

% Import data
results = dlmread('../data/results.txt');
nMC = 1000;
size = 2;

nTimes = length(results) / (size * (nMC + 1));
u = zeros(size, nTimes);

% Time integration parameters.
h = 0.02;
T = 20;
time = 0 : h : T;


labels = {'x', 'y', 'z'};
% Plot deterministic solution
for k = 1 : size
    figure
    hold on
    
    for i = 1 : nTimes
        u(:, i) = results((i - 1) * size + 1 : (i - 1) * size + size);
    end
    plot(time(1:1:end), u(k, 1:1:end), 'r', 'LineWidth', size)
    
    % Plot probabilistic solutions
    for j = size : nMC + 1
        for i = 1 : nTimes
            u(:, i) = results((j - 1) * (size * nTimes) + (i - 1) * size + 1 : ...
                (j - 1) * (size * nTimes) + (i - 1) * size + size);
        end
        plot(time(1:1:end), u(k, 1:1:end), 'b');
    end
    xlabel('t')
    ylabel(labels{k})
end