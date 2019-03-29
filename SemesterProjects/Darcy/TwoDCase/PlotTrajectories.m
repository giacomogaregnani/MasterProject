function PlotTrajectories(M,N,X0,f,g,Time)
% Function used for PLOTTING trajectories

p(1) = plot([-1,-1],[-1,1],'k','LineWidth',2);
hold on
plot([1,1],[-1,1],'k',[-1,1],[-1,-1],'k',[-1,1],[1,1],'k','LineWidth',2)

axis([-1.2 1.2 -1.2 1.2])
axis equal
hold on

W = BrownianMotion2D(Time,N(end),M);

for i = 1:length(N)
    NaiveForPlots(X0,f,g,W(:,1:N(end)/N(i):end),Time,p)
end

% legend(p,{'Domain','big timestep','small timestep'})
end