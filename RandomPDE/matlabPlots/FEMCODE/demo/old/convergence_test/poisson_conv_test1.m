clear mesh sol femspace father;
numIter=10;
offset=2;
mesh= cell(numIter,1);
sol=cell(numIter,1);
femspace = cell(numIter,1);
father = cell(numIter,1);

%% Problem definition
mesh{1} = structured_mesh([0,1,0,1],2);

clear vp;
vp.f = @f1;
vp.elemtype = 'p1';
vp.a = [2, 1; 0.2, 1];

%% Solutions
for i=1:numIter
	[sol{i}, femspace{i}] = poisson(mesh{i}, vp);
    if i<numIter, 
        [mesh{i+1}, father{i+1}] = bisect(mesh{i},'all'); 
    end
end

errL2 = zeros(numIter-1,1);
errH1 = zeros(numIter-1,1);
Ndof = zeros(numIter-1,1);

fat = father{numIter};
for i=numIter-1:-1:1
	errL2(i) = get_L2error(mesh{i}, femspace{i}, sol{i}, ...
		mesh{numIter}, femspace{numIter}, sol{numIter}, fat);
	errH1(i) = get_H1error(mesh{i}, femspace{i}, sol{i}, ...
		mesh{numIter}, femspace{numIter}, sol{numIter}, fat);
	Ndof(i) = femspace{i}.ndof;
	if i>1, fat = father{i}(fat); end
end

%% Convergence rates
figure;
loglog(Ndof(1:end-offset),errL2(1:end-offset),'-x', ...
	[10^2,10^3],[10^(-3), 10^(-4)]);
hold on;
params={'Interpreter','LaTex','FontSize',14};
legend({'$\| u-u_h\|_{L^2}$','$N_{dof}^{-1}$'},params{:});
xlabel('$N_{dof}$',params{:})
ylabel('error in $L^2$-norm',params{:})


figure;
loglog(Ndof(1:end-offset),errH1(1:end-offset),'-x', ...
	[10^2,10^3],[10^(-1.5), 10^(-2)]);
hold on;
legend({'error in H1 norm','slope N^(-1/2)'})
params={'Interpreter','LaTex','FontSize',14};
legend({'$\| u-u_h\|_{H^1}$','$N_{dof}^{-1/2}$'},params{:});
xlabel('$N_{dof}$',params{:})
ylabel('error in $H^1$-norm',params{:})