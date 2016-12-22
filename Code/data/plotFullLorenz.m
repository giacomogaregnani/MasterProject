% Plot trajectory of Lorenz system
results = dlmread('fullResultsLorenz.txt');
t = results(:, 1);
nT = length(t);
Ucom = results(: , 2 : end);

[ha, ~] = tight_subplot(3, 1, [.04 .02], [.1 .05], [.1 .05]);

for i = 1 : 3
    axes(ha(i));
    plot(t, Ucom(:, i));
    if i < 3
        set(gca, 'xticklabel', '')
    end
end
xlabel('t')
