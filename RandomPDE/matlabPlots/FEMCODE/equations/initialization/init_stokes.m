function [vp, mesh] = init_stokes(mesh, femspace, vp)
%INIT_VP Summary of this function goes here
%   Detailed explanation goes here

global global_options;

%% CONSTANTS
dim = size(mesh.elem,2) - 1;

	function initfield(fieldname, value)
	if ~isfield(vp,fieldname)
		vp.(fieldname) = value;
	end
	end

%% QUADRATIC FORMULA FOR STIFFNESS MATRIX (DIFFUSION)
initfield('aquad',max(1,2*(femspace.deg - 1)));
if ~isfield(vp,'alambda') || ~isfield(vp,'aweight')
	[vp.alambda, vp.aweight] = quadpts(dim, vp.aquad);
else
	vp.aquad = NaN; % formula is given, we don't know the deree
end

%% BOUNDARY CONDITIONS
% STRING-SPECIFIED BOUNDARY CONDITIONS CONVERSION
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
    vp.bc = {zeros(1,dim),zeros(1,dim)};
elseif iscell(vp.bc)
    if numel(vp.bc) == 1
        vp.bc = {vp.bc{1}; vp.bc{1}};
    end
elseif isnumeric(vp.bc) 
    vp.bc = {vp.bc; vp.bc};
end

%% PARALLEL COMPUTING (WHERE IMPLEMENTED)
initfield('parallel',false);

%% LINEAR ALGEBRA SOLVER
initfield('solver','backslash');
if strcmp(vp.solver,'iuzawapcg') && ... 
		~(isfield(global_options,'agmg') && global_options.agmg)
	fprintf('AGMG preconditionning unavailable, using backslash');
	vp.solver = 'backslash';
end	

%% VISCOSITY PROPERTIES
if isfield(vp,'a') 
    if ~isnumeric(vp.a) 
        error('The viscosity has to be numeric');
    elseif ~(numel(vp.a) == 1)
        error('The viscosity has to be a single number');
    elseif ~(vp.a > 0)
        error('The viscosity has to be positive')
    end
else
	error('Specify the tensor.')
end
vp.fully_discrete = false; % necessary for assemble_stima to run
vp.symmetric = true;       %


%% QUADRATIC FORMULA FOR RHS
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

%% RHS TYPE (H^-1 not implemented)
[~, vp.NF] = get_force_dim(mesh, vp.f);
vp.fdim = dim;
vp.derivatives = false;        % necessary for assemble_rhs to run
vp.fully_discrete_rhs = false; %
end

