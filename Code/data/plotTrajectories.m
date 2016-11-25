% Plot trajectories (look invariant measures)

nTraj = 10;
y = dlmread('trajHires.txt');
nTimesPerTraj = size(y, 1) / nTraj;
T = 200;
t = 0 : T / nTimesPerTraj : T;

initialCond = [1, zeros(1, 6), 0.0057];

figure 
hold on
for i = 1 : nTraj - 1
    index = (i-1)*nTimesPerTraj + 1;
    plot(t, y(index:index+nTimesPerTraj,5:6));
end