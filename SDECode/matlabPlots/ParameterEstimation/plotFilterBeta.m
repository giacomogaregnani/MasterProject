clc; clear; close all

d = 0.01;
tVec = linspace(0, 1, 1000);

figure
hold on
gray = linspace(0.1, 0.9, 20);
i = 0;
for b = 5 %:length(gray)
    i = i + 1;
    Cb = 1/gamma((b+1)/b);
    k = @(t) Cb / d * exp(-(t / d).^b);
    plot(tVec, k(tVec), 'color', gray(i) * ones(1, 3))
end
ylim = get(gca, 'yLim')

figure
hold on
i = 0;
for b = 50 % :length(gray)
    i = i + 1;
    Cb = 1/gamma((b+1)/b);
    k = @(t) Cb / d^(1/b) * exp(-t.^b / d);
    plot(tVec, k(tVec), 'color', gray(i) * ones(1, 3))
end
% set(gca, 'ylim', ylim)

% d = 1.0;
% figure
% hold on
% i = 0;
% for b = 1:length(gray)
%     i = i + 1;
%     Cb = 1/gamma((b+1)/b);
%     k = @(t) Cb / d^(1/b) * exp(-t.^b / d);
%     plot(tVec, k(tVec), 'color', gray(i) * ones(1, 3))
% end
% set(gca, 'ylim', ylim)