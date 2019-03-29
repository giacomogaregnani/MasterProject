function [vp, mesh] = init_elasticity_orth(mesh, femspace, vp)

dim = size(mesh.elem, 2) - 1;
NT = size(mesh.elem, 1);

function initfield(fieldname, value)
if ~isfield(vp, fieldname)
    vp.(fieldname) = value;
end
end

if ~isfield(vp, 'alambda') || ~isfield(vp, 'aweight')
    initfield('aquad', max(1,2*(femspace.deg - 1)));
    [vp.alambda, vp.aweight] = quadpts(dim, vp.aquad);
else
    initfield('aquad', NaN);
end

if ~isfield(mesh, 'bdflag')
    mesh.bdflag = 'dirichlet';
end

if ~iscell(mesh.bdflag)
    str = mesh.bdflag;
    mesh.bdflag = cell(2,1);
    allfun = @(x) (true(size(x,1), 1));
    allflag = get_bdflag(mesh, allfun);
    if strcmpi(str, 'dirichlet') || strcmpi(str, 'periodic')
        mesh.bdflag{1,1} = allflag;
        mesh.bdflag{2,1} = zeros(0,2);
    elseif strcmpi(str, 'natural') || strcmpi(str, 'neumann')
        mesh.bdflag{1,1} = zeros(0,2);
        mesh.bdflag{2,1} = allflag;
    end
end

if ~isfield(vp, 'bc')
    vp.bc = {zeros(1,dim); zeros(1,dim)};
elseif iscell(vp.bc)
    if numel(vp.bc) == 1
        vp.bc = {vp.bc{1}; vp.bc{1}};
    end
elseif isnumeric(vp.bc)
    vp.bc = {vp.bc; vp.bc};
end

initfield('symmetric',false);
initfield('fully_discrete',false)
NQ = numel(vp.aweight);

initfield('solver', 'backslash');

initfield('fquad',max(1,2*(femspace.deg - 1)));
if vp.fully_discrete
	vp.fquad = vp.aquad;
    vp.flambda = vp.alambda;
    vp.fweight = vp.aweight;
elseif ~isfield(vp,'flambda') || ~isfield(vp,'flambda')
	[vp.flambda, vp.fweight] = quadpts(dim, vp.fquad);
else
	vp.fquad = NaN;
end

if isfield(vp, 'f')
    [~, vp.NF] = get_force_dim(mesh, vp.f);
    vp.fdim = dim;
    if isfield(vp, 'type') && strcmp(vp.type, 'micro_linear_elasticity_orth')
        vp.derivatives = true;
    else
        vp.derivatives = false;
    end
    vp.fully_discrete_rhs = false;
end

end
