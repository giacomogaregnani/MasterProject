
results = dlmread('fullResultsBruss.txt');
t = results(:, 1);
nT = length(t);
Ucom = results(: , 2 : end);
nSpecies = size(Ucom, 2) / 2;
U = Ucom(: , 1 : nSpecies);
V = Ucom(:, nSpecies + 1 : end);
U = [ones(nT, 1), U, ones(nT, 1)];
V = [3 * ones(nT, 1), V, 3 * ones(nT, 1)];

hGrid = 1 / (size(U, 2) - 1);
X = 0 : hGrid : 1;

[XX, TT] = meshgrid(X, t);
surf(XX, TT, U)
xlabel('x')
ylabel('t')