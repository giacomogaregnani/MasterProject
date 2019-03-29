function [sol, femspace, mesh, vp] = poisson_ad(mesh, vp, options)
%POISSON_AD Summary of this function goes here
%   Detailed explanation goes here

%% VP.ADAPT PARSING
if ~isfield(vp,'adapt'), vp.adapt = struct; end
if ~isfield(vp.adapt,'maxit'),   vp.adapt.maxit = 50; end
if ~isfield(vp.adapt,'maxdof'),  vp.adapt.maxdof = 10^6; end
if ~isfield(vp.adapt,'minrelerr'),  vp.adapt.minrelerr = 0; end
if ~isfield(vp.adapt,'minabserr'),  vp.adapt.minabserr = 0; end
if ~isfield(vp.adapt,'marking'), vp.adapt.marking = struct; end

%% OPTIONS PARSING
if nargin<3, options = struct; end
if ~isfield(options,'verbose'), options.verbose = false; end

%% ITERATION
k=1;
while true
    %% SOLVE
	[sol, femspace, mesh, vp] = poisson(mesh, vp); 
    ndof = femspace.ndof; 
    
    %% ESTIMATE
	res = poisson_residual(mesh, vp, femspace, sol);
	abserr = sqrt(sum(res,1));
	relerr = abserr ./ get_H1seminorm(mesh, femspace, sol);
    
    %% DISPLAY   
	if options.verbose
		fprintf('ITER: %d, DOF: %d, REL. ERR: %8.6f, ABS. ERR: %8.6f\n',...
			k, ndof, mean(relerr), mean(abserr));
	end
	
	%% STOP?
	if (k >= vp.adapt.maxit) || (ndof >= vp.adapt.maxdof) || ...
			all(abserr < vp.adapt.minabserr) || ...
			all(relerr < vp.adapt.minrelerr)
		if options.verbose
			fprintf(1,'Iterative process finished.\n');
		end
		break;
	end
	
	%% MARK
	mElements = markelem(res, vp.adapt.marking);
	
	%% REFINE
	mesh = bisect(mesh, mElements);
	k=k+1;
end
end

