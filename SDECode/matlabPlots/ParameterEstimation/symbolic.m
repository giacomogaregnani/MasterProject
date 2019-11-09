clc; clear all; close all
%%

syms t s d b p T
k = 1/(d^(1/b)) * exp(-(t-s)^b/d);
assume(t > 0)
assume(d > 0)
assume(b == 1)

% I = int(k, s, 0, t);
% pretty(simplify(I))
% limit(I, t, inf)

% phi = 1 - I;
% I2 = int(phi, t, 0, T);
% pretty(simplify(I2))
% limit(I2, d, 0)

syms p
assume(p == 1/2)
% assume(p == 2)

f = k * (t-s)^(p);
f = simplify(f);

I = int(f, s, 0, t);
pretty(simplify(I))
limit(I, t, inf)

% syms T
% assume(s == 0)
% f = k * t^p;
% f = simplify(f);
% I = int(f, t, 0, T);
% pretty(simplify(I))