clc; close all; clear;
%% Toda lattice problem n = 3

% ODE definition
f = @(t, y) [-exp(y(4) - y(5)) + exp(y(6) - y(4));   
             exp(y(4) - y(5)) - exp(y(5) - y(6));
             exp(y(5)-y(6)) - exp(y(6) - y(4));
             y(1);
             y(2);
             y(3)];
         
y_0 = [-1.5; 1; 0.5; 1; 2; -1];

[t, y] = ode23s(f, [0, 100], y_0);