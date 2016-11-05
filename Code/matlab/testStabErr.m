results = dlmread('../data/stabErrors.txt');
h = [1 0.5 0.1 0.05 0.01 0.005 0.001 0.0005 0.0001];

loglog(h, results,'o-')
hold on
loglog(h, h, '-k')