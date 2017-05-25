clc; clear; close all

fileNames = {'plotFEM01.txt', 'plotFEM005.txt', 'plotFEM002.txt', 'plotFEM001.txt'};
hTest = [0.1, 0.05, 0.02, 0.01];

for i = 1 : length(fileNames)

    % Read data (exe: PlotPoisson)

    U = dlmread(fileNames{i});
    h = U(1, 1);
    U = U(2:end-1, :);
    UMean = U(end, :);
    
    % Plot
    
    X = 0 : h : 1;
    
    lg = 0.7;
    lightgrey = lg * ones(1, 3);
    
    figure
    plot(X, U, 'color', lightgrey);
    hold on
    p1 = plot(X, UMean, 'r--', 'LineWidth', 2);
    
    % Exact solution
    xEx = linspace(0, 1, 1e5);
    p2 = plot(xEx, 1 / pi^2 * sin(pi * xEx), 'k--', 'LineWidth', 2);
    % p2 = plot(xEx, 1 / pi^2 * (cos(pi * xEx) + 2 * xEx - 1), 'k--', 'LineWidth', 2);
    
    legend([p1, p2], {'mean solution', 'exact solution'})
    title(['h = ', num2str(hTest(i))]);
end