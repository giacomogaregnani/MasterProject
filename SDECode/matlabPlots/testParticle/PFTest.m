clc; clear; close all
%%

obs = dlmread('testSol.txt');
PF = dlmread('test.txt');

t = linspace(0, 1, length(obs));

plot(t, obs, t, PF)