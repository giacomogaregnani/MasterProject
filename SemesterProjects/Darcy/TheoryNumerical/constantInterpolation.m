function [Ux, Uy] = constantInterpolation(f, deltaU)
% Constant interpolation

s = 2 / deltaU;

Ux = zeros(s, s);
Uy = Ux;

for i = 1 : s
    for j = 1 : s
        U = f(-1 + (i+0.5)*deltaU ,-1 + (j+0.5)*deltaU);
        Ux(i, j) = U(1);
        Uy(i, j) = U(2);
    end
end