% Test if the variance depends on h
clc; clear; close all;

MorH = 'h';

results = dlmread('varMCResults_nonZVar.txt');

hEM = results(:, 1);
MEM = results(:, 2);
hUnique = unique(hEM);
hUnique = hUnique(end:-1:1);
MUnique = unique(MEM);
varEM = results(:, 3);
errEM = results(:, 4);

hMin = min(hEM);
MMax = max(MEM);

varTest = varEM(find(hEM == hMin));
errTest = errEM(find(MEM == MMax));

loglog(hUnique, errTest)
hold on
loglog(MUnique.^-1, varTest)
