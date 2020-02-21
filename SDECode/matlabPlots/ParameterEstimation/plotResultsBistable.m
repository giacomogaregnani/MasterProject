clc; clear; close all
addpath('resultsTwo')
%%

hom = dlmread('resultsTwo/ResultsHom_s0.7.txt');
filter = dlmread('resultsTwo/ResultsFilter2_s0.7_b10.txt');
sub = dlmread('resultsTwo/ResultsSub2_s0.7.txt');

aHom = hom(1:2);

zetas = filter(:, 3);

for i = 1 : 2
    
    figure
    hold on
    plot([zetas(1), zetas(end)], [aHom(i), aHom(i)], 'k--');
    plot(zetas, sub(:, i), 'ko-')
    title(['subsampl - parameter ' num2str(i)])
    ylim([0, 1])
    
    figure
    hold on
    plot([zetas(1), zetas(end)], [aHom(i), aHom(i)], 'k--');
    plot(zetas, filter(:, i), 'ko-')
    title(['filter - parameter ' num2str(i)])
    ylim([0, 1])

end