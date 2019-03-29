clc; clear; close all


syms h H p

f = (H^2 - h^2)^2;
% f = (H^);
assume(h > 0);
assume(H > 0);

I = int(f/(2*h^p), H, h-h^p, h+h^p)
