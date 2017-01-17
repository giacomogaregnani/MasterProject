clc; clear; close all

h = 0.01;
hText = floor(h * 1e5);
nExp = length(h);

phi = @(Y) sum(Y.^2, 2);

YMP = dlmread(['invMeasMPLorenz', num2str(hText(1)), '.txt']);
phiMP = phi(YMP);

scatter3(YMP(1:100, 1), YMP(1:100, 2), YMP(1:100, 3), '.')

yFull = dlmread(['invMeasMPLorenzFull', num2str(hText(1)), '.txt']);
M = 100;


lightgray = [0.8, 0.8, 0.8]
index = 1;
selection = 1 : 5 : 4001;

figure
hold on
for i = 1 : 20
    y = yFull(index : index + 4000, :);
    t = y(:, 1);
    y = y(:, 2:4);
    index = index + 4001;
    if i == 20
        plot(t(selection), y(selection, 1), 'k', 'LineWidth', 2);
    else
        plot(t(selection), y(selection, 1), 'color', lightgray);
    end
end
xlim([0, 40.5])
xlabel('t')