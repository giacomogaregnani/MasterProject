clear; close all; clc;

Y = dlmread('fullSolution.txt');
N = 100;

figure
hold on

index = 1;
times = linspace(0.1, 10, N);

lightgray = [0.8, 0.8, 0.8];
for i = 1 : N
    solution = Y(index : index + N - 1, 2 : end); 
    index = index + N;
    
    plot(times, solution, 'color', lightgray)
end

index = 1;
for i = 1 : N
    times = Y(index : index + N - 1, 1);
    solution = Y(index : index + N - 1, 2 : end); 
    index = index + N;
    
    plot(times, solution, 'k', 'LineWidth', 2)
end
