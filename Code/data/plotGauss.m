% Plot trajectories (look invariant measures)

clear
y = dlmread('meansHires.txt');
V = dlmread('varHires.txt');
sizeODE = 2;
% Plot
T = 50;
h = 0.001;
t = 0 : h : T;
nData = length(t);

i = 1;
varEl = zeros(nData, sizeODE);
for j = 1 : nData 
   var = V(i:i+sizeODE-1, :);
   varEl(j, :) = 2.0 * sqrt(abs(diag(var)));
   i = i + sizeODE - 1;
end

plot(t, y, 'LineWidth', 2);
hold on
plot(t, y + varEl, 'r');
