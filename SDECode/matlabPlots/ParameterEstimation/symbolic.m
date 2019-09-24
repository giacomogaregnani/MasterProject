clc; clear all; close all
%%

syms t s a d b c
k = 1/(d^(1/b)) * exp(-(t-s)^b/d);
dtk = (t-s)^a*diff(k, t);
dtk = simplify(dtk)

assume(t, 'real')
assume(d, 'real')
assume(d > 0)
assume(a == 3/2)
assume(b == 1)
% assume(c > 0)
% assume(c == 1)

I = int(dtk, s, 0, t)

limit(I, t, inf)