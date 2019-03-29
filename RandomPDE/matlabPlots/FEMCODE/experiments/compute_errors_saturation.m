function compute_errors_saturation(experiment,p,q)
% ERROR COMPUTATION FOR FE-HMM PROBLEMS

prefix = 'iter';


	

direc = ['~/repos/experiments/comp/' experiment '/'];
rhs =  str2func(['rhs_' experiment]);
ahom = str2func(['ahom_' experiment]);

if q==1
	umicro = 'p1b';
else
	umicro = ['p' num2str(q)];
end
pmicro = ['p' num2str(max(1,q-1))];
macro =  ['p' num2str(p)];
exp_name = [macro 'vs' umicro pmicro 'sat800']; % SATURATION !
direcout = [direc exp_name '/']; 
if ~exist(direcout,'dir'), mkdir(direcout); end


finished = 0;
for i=1:200
	if exist([direcout prefix num2str(i) '.mat'],'file')
		finished = i;
	end
end

load([direcout prefix num2str(finished) '.mat'], ...
	'macMesh','macVp','macSol','macFemspace');

% FINE SOLUTION
vpf.elemtype = ['p' num2str(min(p+1,3))];
vpf.a = ahom;
vpf.f = rhs;
vpf.bc = {0,0};
vpf.solver = 'agmg';
vpf.symmetric = true;


% INITIALIZATION
ndof   = zeros(finished,1);
ndofcomp = zeros(finished,1);
h1semierr = zeros(finished,1);
elemh1semierr = cell(finished,1);



for i=1:finished 
  
  load([direcout prefix num2str(i) '.mat'], ...
	'macMesh','macVp','macSol','macFemspace');

	fprintf('Processing %d/%d : refine,',i,finished);

	% Two uniform refinements
	[meshf,  fat1] = uniformrefine(macMesh);
	[meshf,  fat2] = uniformrefine(meshf);
	fat = fat1(fat2); 
	
	% fine solution
	fprintf('F, ');
	[solf, femspacef, meshf, vpf] = poisson(meshf, vpf);
	
	%% error computation
	NT = size(macMesh.elem, 1);
	
	fprintf('E, ');
	[ ~, ~, h1semierr(i), ~, ~, elemh1semierr{i}  ] = ...
		get_H1error(macMesh, macFemspace, macSol, ...
		meshf, femspacef, solf, fat);
  
  ndof(i) = macFemspace.ndof;
  ndofcomp(i) = ndof(i) + sum(macVp.micNdof(:));
  
	% FINISHED
    fprintf('\n');
end


[H1norm, L2norm, H1seminorm] = get_H1norm(meshf, femspacef, solf);

%% save error computation
save([direcout prefix 'error' num2str(i) '.mat'], ...
	'ndof','ndofcomp', ...
	'h1semierr','H1norm','L2norm','H1seminorm');
end
