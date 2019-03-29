function [Ux, Uy] = linearInterpolation(f, deltaU)
% Linear interpolation 

s = 2 / deltaU;

Ux = zeros(s + 1, s + 1);
Uy = Ux;

for i = 1 : s + 1
    for j = 1 : s + 1
        U = f(-1 + (i-1)*deltaU ,-1 + (j-1)*deltaU);
        Ux(i, j) = U(1);
        Uy(i, j) = U(2);
    end
end