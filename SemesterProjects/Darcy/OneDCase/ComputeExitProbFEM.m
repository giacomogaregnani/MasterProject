function phi = ComputeExitProbFEM(X0)
% CMCS : Compute EXIT PROBABILITY using 1D FEM.

%  Copyright (C) EPFL, CMCS

addpath('CMCS_FEM1D')

a            = 1;
N            = 6/0.005;
dt           = 0.01;

fem = 'P1';
Theta = 0.5;

for i = 1 : length(dt)    
    nx = N(i);      
    [vertices,elements,boundaries] = mesh1D(-a, a , nx);

    solution = Parabolic1D_Solver(elements,vertices,boundaries,fem,'FEM_DATA',[dt(i) Theta]);
end

plot(linspace(-1, 1, N +1), solution(:, end))

phi = interp1(linspace(-1, 1, N + 1), solution(:, end), X0);



