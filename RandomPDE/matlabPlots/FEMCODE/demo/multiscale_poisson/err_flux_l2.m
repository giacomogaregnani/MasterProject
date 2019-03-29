function out = err_flux_l2(f1,f2)
N1 = size(f1,1);
N2 = size(f2,1);
x1 = linspace(0,1,N1);
x2 = linspace(0,1,N2);
f2 = interp1(x2, f2, x1)';
out = sqrt(trapz(x1, (f1-f2).^2));