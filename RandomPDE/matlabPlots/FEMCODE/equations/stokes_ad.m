function [usol, ufemspace, psol, pfemspace, mesh, vp] = ...
	stokes_ad(mesh, vp, options)
%STOKES_AD solves adaptively a Stokes problem given by mesh and vp
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

%% CONSTANTS
dim = size(mesh.elem,2)-1;

%% ITERATION
k=1;
while true
	%% SOLVE
	[usol, ufemspace, psol, pfemspace, mesh, vp] = ...
		stokes(mesh, vp, options);
    ndof = dim*ufemspace.ndof + pfemspace.ndof;
	
	%% ESTIMATE
	res = stokes_residual(mesh, vp, ufemspace, usol, pfemspace, psol);
	abserr = sqrt(sum(res,1));
	relerr = abserr ./ get_H1seminorm(mesh, ufemspace, usol);
	
	%% DISPLAY
	if options.verbose
		fprintf('ITER: %d, DOF: %d, REL. ERR: %8.6f, ABS. ERR: %8.6f\n', ...
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

