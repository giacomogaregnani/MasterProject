function ErrorPlots(errBernoulliPhi,errNaivePhi,errBernoulli,errNaive,N,Time)
% Plot the error for orders analysis

h = (Time(2)-Time(1))./N;
IndForPlots = ceil(length(h)/2) + 1;
figure
loglog(h,errNaive,'ro-')
hold on
loglog(h,errBernoulli,'b*-')
loglog(h,sqrt(h)*(errNaive(IndForPlots)/sqrt(h(IndForPlots))),'k--')
loglog(h,h*(errBernoulli(IndForPlots)/h(IndForPlots)),'k')
grid on
h_legend = legend('DEM','CEM','h^{0.5}','h');
set(h_legend,'Location','best','FontSize',13);
xlabel('h')
title('Convergence of the exit time')

figure
loglog(h,errNaivePhi,'ro-')
hold on
loglog(h,errBernoulliPhi,'b*-')
loglog(h,sqrt(h)*(errNaivePhi(IndForPlots)/sqrt(h(IndForPlots))),'k--')
loglog(h,h*(errBernoulliPhi(IndForPlots)/h(IndForPlots)),'k')
grid on
h_legend = legend('DEM','CEM','h^{0.5}','h');
set(h_legend,'Location','best','FontSize',13);
xlabel('h')
title('Convergence of the exit probability')

end

