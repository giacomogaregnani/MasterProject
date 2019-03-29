function [vp, mesh] = init_elliptic(mesh, femspace, vp)
%INIT_VP Summary of this function goes here
%   Detailed explanation goes here

global global_options;

%% CONSTANTS
dim = size(mesh.elem,2) - 1;
NT = size(mesh.elem,1);
LD = femspace.ldof;

	function initfield(fieldname, value)
	if ~isfield(vp,fieldname)
		vp.(fieldname) = value;
	end
	end

%% QUADRATIC FORMULA FOR STIFFNESS MATRIX (DIFFUSION)
if ~isfield(vp,'alambda') || ~isfield(vp,'aweight')
    initfield('aquad',max(1,2*(femspace.deg - 1)));
    [vp.alambda, vp.aweight] = quadpts(dim, vp.aquad);
else
	initfield('aquad',NaN);
end

%% QUADRATIC FORMULA FOR ADVECTION TERMS
if isfield(vp,'b')
	initfield('bquad',max(1,2*femspace.deg - 1));
	if ~isfield(vp,'blambda') || ~isfield(vp,'bweight')
		[vp.blambda, vp.bweight] = quadpts(dim, vp.bquad);
	else
		vp.bquad = NaN;  % formula is given, we don't know the deree
	end
end

%% QUADRATIC FORMULA FOR REACTION TERMS
if isfield(vp,'c')
	initfield('cquad',max(1,2*femspace.deg));
	if ~isfield(vp,'clambda') || ~isfield(vp,'cweight')
		[vp.clambda, vp.cweight] = quadpts(dim, vp.cquad);
	else
		vp.cquad = NaN;  % formula is given, we don't know the deree
	end
end

%% BOUNDARY CONDITIONS
% STRING-SPECIFIED BOUNDARY FLAG CONVERSION
if ~isfield(mesh,'bdflag')
    mesh.bdflag = 'dirichlet';
    fprintf('Boundary type not specified, using Dirichlet\n');
end
if ~iscell(mesh.bdflag)
    str = mesh.bdflag;
    mesh.bdflag = cell(2,1);
	allfun = @(x)(true(size(x,1),1));
	allflag = get_bdflag(mesh,allfun);
	if strcmpi(str, 'dirichlet') || strcmpi(str, 'periodic')    
		mesh.bdflag{1} = allflag;
		mesh.bdflag{2} = zeros(0,2);
	elseif strcmpi(str, 'natural') || strcmpi(str, 'neumann')
		mesh.bdflag{1} = zeros(0,2);
		mesh.bdflag{2} = allflag;
	end
end
if ~isfield(vp,'bc')
    vp.bc = num2cell(zeros(2,1));
elseif iscell(vp.bc)
    if numel(vp.bc) == 1
        vp.bc = {vp.bc{1}; vp.bc{1}};
    end
elseif isnumeric(vp.bc) && (numel(vp.bc) == 2)
    vp.bc = {vp.bc(1); vp.bc(2)};
else
    vp.bc = {vp.bc; vp.bc};
end

%% PARALLEL COMPUTING (WHERE IMPLEMENTED)
initfield('parallel',false);

%% LINEAR ALGEBRA SOLVER
initfield('solver','backslash');
if strcmp(vp.solver,'agmg') && ... 
		~(isfield(global_options,'agmg') && global_options.agmg)
	fprintf('AGMG solver unavailable, using backslash');
	vp.solver = 'backslash';
end	

%% TENSOR PROPERTIES
initfield('symmetric',false);
initfield('fully_discrete',false)
NQ = numel(vp.aweight);
if isfield(vp,'a') && isnumeric(vp.a)
	if (numel(vp.a) == 1)
		if vp.a <= 0
			error('tensor you specified is not positive definite');
		end
		vp.a=eye(dim)*vp.a;
		vp.symmetric = true;
	elseif (numel(vp.a) == dim)
		if any(vp.a <= 0)
			error('tensor you specified is not positive definite');
		end
		vp.a=diag(vp.a(:));
		vp.symmetric =true;
	elseif (numel(vp.a) == dim^2)
		if (sum(sum(abs(vp.a-vp.a'))) == 0)
			vp.symmetric =true;
		end
    else
        if size(vp.a,1) == NT && size(vp.a,2) == dim && ...
                size(vp.a,3)== dim && size(vp.a,4) == NQ
            vp.fully_discrete = true;
        else
            error('Unknown type of fully discrete data');
        end
	end
elseif ~isfield(vp,'a')
	error('Specify the tensor.');
end

% %% ADVECTION PROPERTIES
% if ~isfield(vp,'b')
%     display('VELOCITY FIELD NOT SPECIFIED, USING zeros(dim,1) instead.');
% end
% 
% %% REACTION PROPERTIES
% if ~isfield(vp,'c')
%     display('REACTION TERM NOT SPECIFIED, USING 0 instead.');
% end

%% QUADRATIC FORMULA FOR RHS
%
% TODO: ADAPTIVELY CHOOSE FOR ADVECTION, REACTION DOMINATED PROBLEMS
%
if (isfield(vp,'b') || isfield(vp,'c')) && ~isfield(vp,'fquad');
    display('Quadrature formula for RHS might be not appropriate, if not diffusion dominated problem.');
end

initfield('fquad',max(1,2*(femspace.deg - 1)));
if vp.fully_discrete
	vp.fquad = vp.aquad;
    vp.flambda = vp.alambda;
    vp.fweight = vp.aweight;
elseif ~isfield(vp,'flambda') || ~isfield(vp,'flambda')
	[vp.flambda, vp.fweight] = quadpts(dim, vp.fquad);
else
	vp.fquad = NaN; % quad. formula was given, order not known
end

%% RHS TYPE (DERIVATIVES OR NOT)
vp.fdim = 1;
[FD, vp.NF] = get_force_dim(mesh, vp.f);
if (FD==1)
	vp.derivatives=false; 
elseif (FD == 1 + dim)
	vp.derivatives=true;
else
	error('wrong size of the Force field');
end
if isnumeric(vp.f) 
    if size(vp.f,1) == 1 && size(vp.f,2) == FD && size(vp.f,3) == vp.NF
        vp.fully_discrete_rhs = false;
    elseif (size(vp.f,1) == NT) && (size(vp.f,4) == NQ)
        vp.fully_discrete_rhs = true;
    else
        error('Unknown numeric force description');
    end
end
end

