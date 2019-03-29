%% generate mesh
meshs = structured_mesh([0,1,0,1],[3,3],struct('bc','periodic'));
meshs.elem([9,10],:) = [];
meshs = uniformrefine(meshs);
dim = size(meshs.elem,2) - 1;

%% set variational problem
clear vp;
vp.f = @fstokesmicro;
vp.elemtype = 'p2';
vp.pelemtype = 'p1';
vp.a = 1;

%% Solve
numIterTotal=10; % total number of iterations
numIterBisection = 4; % number of bisections in the end

%% initialization
[~, NF] = get_force_dim(meshs,vp.f);
usol = cell(numIterTotal,1);
psol = cell(numIterTotal,1);
mesh = cell(numIterTotal,1); 
allRresiduals = cell(numIterTotal,1); 
mesh{1} = meshs;
ufemspace = cell(numIterTotal,1);
pfemspace = cell(numIterTotal,1);
father = cell(numIterTotal,1);
errL2 = zeros(numIterTotal-1,1,NF);
errH1 = zeros(numIterTotal-1,1,NF);
errRes = zeros(numIterTotal-1,1,NF);
Ndof  = zeros(numIterTotal-1,1);
allResiduals = cell(numIterTotal-numIterBisection-1,1);

%% ADAPTIVE REFINEMENT PHASE
for i=1:numIterTotal-numIterBisection-1
	[usol{i}, ufemspace{i}, psol{i}, pfemspace{i}, mesh{i}] = ...
		stokes(mesh{i}, vp);
	res = stokes_residual(mesh{i}, vp, ufemspace{i}, usol{i}, ...
		pfemspace{i}, psol{i});
	
	errRes(i,:,:) = sqrt(sum(res,1));
	allResiduals{i}=res;
	res = sqrt(res);
	
	markedElements = markelem(res);
	
	fprintf(1,'Adaptive refinement iteration %d / %d\n', ...
		i,numIterTotal-numIterBisection-1);
	[mesh{i+1}, father{i+1}]= bisect(mesh{i}, markedElements);
end

%% BISECTION REFINEMENT PHASE
for i= numIterTotal-numIterBisection:numIterTotal
	[usol{i}, ufemspace{i}, psol{i}, pfemspace{i}, mesh{i}] = ...
		stokes(mesh{i}, vp);
	if i<numIterTotal
		fprintf(1,'Bisection %d / %d\n', ...
			numIterBisection + i - numIterTotal +1, numIterBisection);
		[mesh{i+1}, father{i+1}]= bisect(mesh{i},'all');
	end
end

%% ERROR COMPUTATION
fat = father{numIterTotal};
for i=numIterTotal-1:-1:1
	fprintf(1,'Error computation %d / %d\n', ...
		numIterTotal-i, numIterTotal-1);
	[errH1(i,:,:), errL2(i,:,:)] = ...
		get_H1error(mesh{i}, ufemspace{i}, usol{i}, ...
		mesh{numIterTotal}, ufemspace{numIterTotal}, ...
		usol{numIterTotal}, fat);
	Ndof(i) = ufemspace{i}.ndof * dim + pfemspace{i}.ndof;
	if i>1, 
		fat = father{i}(fat); 
	end
end

%% Convergence rates in L2 norm
figure;
loglog(Ndof(1:end-numIterBisection+1), ...
	errL2(1:end-numIterBisection+1,:,1),'-x', ...
	Ndof(1:end-numIterBisection+1), ...
	errL2(1:end-numIterBisection+1,:,2),'-x', ...
	Ndof(1:end-numIterBisection+1), ...
	errRes(1:end-numIterBisection+1,:,1).^2/150,'-x', ...
	Ndof(1:end-numIterBisection+1), ...
	errRes(1:end-numIterBisection+1,:,2).^2/150,'-x', ...
	[10^2,10^5],[10^(-3), 10^(-6)]);
hold on;
params={'Interpreter','LaTex','FontSize',12};
legend({'$\| u^1-u^{1,h}\|_{L^2}$', '$\| u^2-u^{2,h}\|_{L^2}$', ...
	'$(e^1_{estimate})^2/150$', '$(e^2_{estimate})^2/150$', ...
	'$N_{dof}^{-1}$'}, ...
	params{:});
xlabel('$N_{dof}$',params{:})
ylabel('errors in $L^2$-norm',params{:})

%% Convergence rates in H1 norm
figure;
loglog(Ndof(1:end-numIterBisection+1),...
	errH1(1:end-numIterBisection+1,:,1), '-x', ...
	Ndof(1:end-numIterBisection+1),...
	errH1(1:end-numIterBisection+1,:,2), '-x', ...
	Ndof(1:end-numIterBisection+1), ...
	errRes(1:end-numIterBisection+1,:,1),'-x', ...
	Ndof(1:end-numIterBisection+1), ...
	errRes(1:end-numIterBisection+1,:,2),'-x', ...
	[10^2,10^5],[10^(-1.5), 10^(-3)]);
hold on;
params={'Interpreter','LaTex','FontSize',12};
legend({'$\| u^1-u^{1,h}\|_{H^1}$','$\| u^1-u^{1,h}\|_{H^1}$', ...
	'$e^1_{estimate}$','$e^2_{estimate}$', ...
	'$N_{dof}^{-1/2}$'}, ...
	params{:});
xlabel('$N_{dof}$',params{:})
ylabel('error in $H^1$-norm',params{:})


%% final mesh plot
figure;
simpplot(mesh{numIterTotal-numIterBisection});

%% plot pressure
figure; 
simpplot_sol(mesh{numIterTotal-numIterBisection},psol{numIterTotal-numIterBisection});

%% plot velocity
figure; 
simpplot_sol(mesh{numIterTotal-numIterBisection},usol{numIterTotal-numIterBisection}(:,1));


