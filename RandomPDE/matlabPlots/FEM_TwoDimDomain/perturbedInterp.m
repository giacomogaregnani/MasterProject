clc; clear; close all

f = @(x) cos(2 * pi * x(1, :))' .* sin(2 * pi * x(2, :))';

NRef = 400;
[vRef, bRef, eRef] = poimesh('squareg', NRef, NRef);

% useless
data   = read_DataFile('InverseData', 2, []);

fRef = f(vRef);

NVec = 2.^[2:8]
err = [];

for N = NVec
    [v, b, e] = poimesh('squareg', N, N);
    [~, ~, mesh]  = Elliptic_Solver(2, e, v, b, 'P1', data, [], [], false);

    
    fN = f(v);
    fEval = evalFEM(mesh, fN, vRef);
    err = [err max(abs(fEval - fRef))];
end

loglog(1./NVec, err, 'ko--')
hold on
loglog(1./NVec, 100./NVec.^2, 'k')