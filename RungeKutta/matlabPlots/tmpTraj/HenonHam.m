clear; clc; close all;

%% Load data
nMC = 10; h = 0.01;
allSolutionsStep = dlmread('../IMChaosStepLong.txt');
totSteps = size(allSolutionsStep, 1);
eachSteps = totSteps / nMC;
allSolutionsAdd = dlmread('../IMChaosAddLong.txt');

%% Plots settings
enhanced = 1;
fontsizeLAB = getLatexTextSize('normalsize', 'enhanced', enhanced);
fontsizeTICK = getLatexTextSize('small', 'enhanced', enhanced);
WTraj = 7; WHor = 9; H = 4; offset = 0.06;

%% Trajectory step

lightgrey = [0.5, 0.5, 0.5];
darkgrey = [0.3, 0.3, 0.3];

indexEnd = 35000;

fig = createFigure(WTraj, H, 'enhanced',enhanced);
plot(allSolutionsStep(1 : indexEnd, 3), ...
    allSolutionsStep(1 : indexEnd, 4), ...
    'color', darkgrey)
set(gca, 'fontsize', fontsizeTICK);
set(get(gca, 'xlabel'), 'fontsize', fontsizeLAB);
set(get(gca, 'ylabel'), 'fontsize', fontsizeLAB);
axis equal
ylim([-0.5, 0.7])
xlim([-0.6, 0.6])
set(gca, 'yticklabel', [])
set(gca, 'xticklabel', [])
ax = axis;

print -depsc2 ../../../Reports/Leysin2017/HHStep.eps


%% Trajectory add

fig = createFigure(WTraj, H, 'enhanced',enhanced);
plot(allSolutionsAdd(1 : indexEnd, 3), ...
    allSolutionsAdd(1 : indexEnd, 4), ...
    'color', darkgrey)
set(gca, 'fontsize', fontsizeTICK);
set(get(gca, 'xlabel'), 'fontsize', fontsizeLAB);
set(get(gca, 'ylabel'), 'fontsize', fontsizeLAB);
axis equal
ylim([-0.5, 0.7])
xlim([-0.6, 0.6])
set(gca, 'yticklabel', [])
set(gca, 'xticklabel', [])
ax = axis;

print -depsc2 ../../../Reports/Leysin2017/HHAdd.eps


%% Chaos on one component

time = linspace(0, h * eachSteps, eachSteps);

fig = createFigure(WHor, H, 'enhanced',enhanced);
hold on
for i = 1 : 3
    thisTrajectory = allSolutionsStep(eachSteps * (i - 1) + 1 : ...
        eachSteps * i, :);
    plot(time, thisTrajectory(:, 3));%, 'color', lightgrey);
end

axis tight

ylabel('$q_1$', 'interpreter', 'LaTex')
box on
set(gca, 'fontsize', fontsizeTICK);
set(get(gca, 'xlabel'), 'fontsize', fontsizeLAB);
set(get(gca, 'ylabel'), 'fontsize', fontsizeLAB);
xlim([400, 600])
xlabel('$t$', 'interpreter', 'LaTex')

print -depsc2 ../../../Reports/Leysin2017/HHStepChaos.eps

xlim([0, 200])
print -depsc2 ../../../Reports/Leysin2017/HHStepChaos2.eps



