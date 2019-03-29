function [xout,yout,stats] = iuzawapcg(A,B,f,M, options)
%%IUZAWAPCG inexact Uzawa preconditioned conjugate gradient method
% 
% This method solves a symmetric system
%      |A  B'| |u|  = |f|
%      |B  0 | |p|  = |0|
% where A is positive definite. 
%
% ALGORITHM: third algorithm from
% Peters, Reichelt, and Reusken : Fast Iterative Solvers for Discrete
% Stokes Equations, SIAM Journal of Scientific Computing, (1-21) 2005
ticoverall = tic;
if nargin<4, error('At least four parameters are necessary'); end
if nargin<5, options=struct; end
if ~isfield(options,'verbose'), options.verbose = false; end
if ~isfield(options,'adaptive'), options.adaptive = false; end
if ~isfield(options,'numInnerIter'), options.numInnerIter = 3; end
if ~isfield(options,'globalTolerance'), options.globalTolerance = 1e-10; end

if options.verbose
	fprintf('### INEXACT UZAWA PCG START\n');
end


%% ALGORITHM SETTINGS AND TOLERANCES
globalTolerance = options.globalTolerance;
precTolerance = 0.5;
maxOuterIter = 200;
numInnerIter = options.numInnerIter;

%% PRECONDITIONING SETUP
prec = full(diag(M)).^(-1);
agmg(A,1,1,[],[],[],[],1);

%% GLOBAL INITIALIZATION
N1 = size(B,2);
N2 = size(B,1);
NF = size(f,3);
xout = zeros(N1, NF);
yout = zeros(N2, NF);
stats.errors = zeros(NF,1);
stats.iters = zeros(NF,1);
stats.prec = zeros(NF,1);
stats.jobtimes = zeros(NF,1);

for j=1:NF
ticjob = tic;

%% INITIALIZATION 
x = zeros(N1,1);
y = zeros(N2,1);
err = 1;
noninc = 0;
r1 = f(:,:,j) - A*x - B'*y;
k = 0;

if options.verbose >= 1
	fprintf('## JOB %d/%d, DOF: %8.0u, INNER ITER: %d\n', ...
		j, NF, N1+N2, numInnerIter);
end


%% THE CORE
while (err>globalTolerance) && (k<maxOuterIter)
	k = k+1;
	if options.verbose >=2
		fprintf('# ITER: %d, (', k);
	end
	xold = x;
	w = x + agmg(A,r1,1,precTolerance,[],[],[],3);

	% SOLVING of z=Phi(B*w) approximately by several steps of PCG
	xx = zeros(N2,1); % solution
	xxbar = zeros(N1,1);
	xxhat = zeros(N1,1);
	bb = B*w;         % right hand side
	%rr = bb - B*(A\(B'*xx));
	[pomrr, ~] = agmg(A,B'*xx,1,precTolerance,[],[],[],3);
	rr = bb - B*pomrr;
	zz = prec .* rr;
	pp = zz;
	for kk=0:numInnerIter -1
		pom1 = B'*pp;
		[pom2, ~] = agmg(A,pom1,1,precTolerance,[],[],[],3);
		pom = B*pom2; % aux
		zzrr = dot(zz,rr);   % aux
		alpha = zzrr/dot(pp,pom);
		xx = xx + alpha*pp;
		xxbar = xxbar + alpha*pom1;
		xxhat = xxhat + alpha*pom2;
		rr = rr - alpha*pom;
		zz = prec .* rr;
		beta = dot(zz,rr)/zzrr;
		pp = zz + beta*pp;
		if options.verbose >= 2, fprintf('-'); end
	end
	if options.verbose >= 2, fprintf('), '); end
	z = xx;
	zhat = xxhat;
	zbar = xxbar;

	x = w - zhat;
	y = y + z;
	r1 = r1 - A*(x-xold) - zbar;
	
	lasterr = err;
    err = sqrt((x-xold)'*A*(x-xold))/sqrt(x'*A*x);
	if options.verbose >= 2
		fprintf('ERR: %12.8g\n', err);
	end
	
	if options.adaptive && (k > 1)
		if err > 0.7 * lasterr
			numInnerIter = numInnerIter + 1;
			noninc = 0;
		else
			noninc = noninc + 1;
		end
		if (numInnerIter > 3) && ...
				((err < 0.2 * lasterr) || (noninc > 2))
			numInnerIter = numInnerIter - 1;
		end
	end
end
if options.verbose >= 1
	fprintf('AGMG PRECONDITIONNER CALLED %d TIMES\n', ...
		numInnerIter*k);
end
xout(:,j) = x;
yout(:,j) = y;
stats.errors(j) = err;
stats.iters(j) = k;
stats.prec(j) = (numInnerIter+1)*k;
stats.jobtimes(j) = toc(ticjob);
end
stats.overall = toc(ticoverall);
