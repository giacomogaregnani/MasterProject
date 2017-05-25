clear; clc; close all;

%% Load data
nMC = 10; h = 0.01;
allSolutionsStep = dlmread('IMChaosStepLong.txt');
totSteps = size(allSolutionsStep, 1);
eachSteps = totSteps / nMC;
allSolutionsAdd = dlmread('IMChaosAddLong.txt');

%% Plots settings
enhanced = 1;
fontsizeLAB = getLatexTextSize('normalsize', 'enhanced', enhanced);
fontsizeTICK = getLatexTextSize('small', 'enhanced', enhanced);
WTraj = 6; WHor = 14; H = 6; offset = 0.05;

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
ylim([-0.7, 0.7])
xlim([-0.6, 0.6])
set(gca, 'xtick', [-0.5, 0.5])
set(gca, 'ytick', [-0.5, 0.5])
set(gca, 'yticklabel', [])
set(gca, 'xticklabel', [])
ax = axis;
text(ax(1) + offset, -0.5, ['-0.5'],...
    'HorizontalAlignment','Left', 'fontname', 'CMUSerif', 'fontsize', fontsizeTICK);
text(ax(1) + offset, 0.5, ['0.5'],...
    'HorizontalAlignment','Left', 'fontname', 'CMUSerif', 'fontsize', fontsizeTICK);
text(-0.5, ax(3) + 2 * offset, -0.5, ['-0.5'],...
    'HorizontalAlignment','Center', 'fontname', 'CMUSerif', 'fontsize', fontsizeTICK);
text(0.5, ax(3) + 2 * offset, 0.5, ['0.5'],...
    'HorizontalAlignment','Center', 'fontname', 'CMUSerif', 'fontsize', fontsizeTICK);

print -depsc2 ../../Reports/RandTimeStep/VERSION3/HHStep.eps


%% Trajectory add

fig = createFigure(WTraj, H, 'enhanced',enhanced);
plot(allSolutionsAdd(1 : indexEnd, 3), ...
    allSolutionsAdd(1 : indexEnd, 4), ...
    'color', darkgrey)
set(gca, 'fontsize', fontsizeTICK);
set(get(gca, 'xlabel'), 'fontsize', fontsizeLAB);
set(get(gca, 'ylabel'), 'fontsize', fontsizeLAB);
axis equal
ylim([-0.7, 0.7])
xlim([-0.6, 0.6])
set(gca, 'xtick', [-0.5, 0.5])
set(gca, 'ytick', [-0.5, 0.5])
set(gca, 'yticklabel', [])
set(gca, 'xticklabel', [])
ax = axis;
text(ax(1) + offset, -0.5, ['-0.5'],...
    'HorizontalAlignment','Left', 'fontname', 'CMUSerif', 'fontsize', fontsizeTICK);
text(ax(1) + offset, 0.5, ['0.5'],...
    'HorizontalAlignment','Left', 'fontname', 'CMUSerif', 'fontsize', fontsizeTICK);
text(-0.5, ax(3) + 2 * offset, -0.5, ['-0.5'],...
    'HorizontalAlignment','Center', 'fontname', 'CMUSerif', 'fontsize', fontsizeTICK);
text(0.5, ax(3) + 2 * offset, 0.5, ['0.5'],...
    'HorizontalAlignment','Center', 'fontname', 'CMUSerif', 'fontsize', fontsizeTICK);

print -depsc2 ../../Reports/RandTimeStep/VERSION3/HHAdd.eps


%% Chaos on one component

time = linspace(0, h * eachSteps, eachSteps);

fig = createFigure(WHor, H, 'enhanced',enhanced);
hold on
for i = 1 : nMC
    thisTrajectory = allSolutionsStep(eachSteps * (i - 1) + 1 : ...
        eachSteps * i, :);
    plot(time, thisTrajectory(:, 3), 'color', lightgrey);
end

axis tight
xlim([0, 600])

ylabel('$q_1$', 'interpreter', 'LaTex')
box on
set(gca, 'fontsize', fontsizeTICK);
set(get(gca, 'xlabel'), 'fontsize', fontsizeLAB);
set(get(gca, 'ylabel'), 'fontsize', fontsizeLAB);
set(gca, 'xtick', [])
set(gca, 'xticklabel', [])

print -depsc2 ../../Reports/RandTimeStep/VERSION3/HHStepChaos.eps


%% Plot the Hamiltonian for each trajectory

fig = createFigure(WHor, H, 'enhanced',enhanced);

hamError = zeros(length(time) - 1, nMC);
for i = 1 : nMC
    hamError = abs(allSolutionsAdd(eachSteps * (i - 1) + 2 : eachSteps * i, 5) - ...
        allSolutionsAdd(eachSteps * (i - 1) + 1, 5));
   h(i) = semilogy(time(2 : end), hamError, 'color', lightgrey);
   hold on     
end

for i = 1 : nMC
    hamError = abs(allSolutionsStep(eachSteps * (i - 1) + 2 : eachSteps * i, 5) - ...
        allSolutionsStep(eachSteps * (i - 1) + 1, 5));
    h(nMC + i) = semilogy(time(2 : end), hamError, 'color', 'k');
end

xlim([0, 600])
xlabel('$t$', 'interpreter', 'LaTex')
ylabel('$|H(p,q) - H(p_0, q_0)|$', 'interpreter', 'LaTex')
legend([h(1), h(nMC + 1)], {'Additive noise', 'Random time step'}, 'Location', 'best')
set(gca, 'fontsize', fontsizeTICK);
set(get(gca, 'xlabel'), 'fontsize', fontsizeLAB);
set(get(gca, 'ylabel'), 'fontsize', fontsizeLAB);
box on

print -depsc2 ../../Reports/RandTimeStep/VERSION3/HHStepHam.eps

