function [sol, femspace, mesh, vp] = poisson_dg(mesh, vp, options)
%POISSON_DG solves a Poisson problem with DG given by mesh and vp

%% OPTIONS PARSING
if nargin<3, options = struct; end
if ~isfield(options,'verbose'), options.verbose = false; end

%% FINITE ELEMENTS INITIALIZATION
mesh = gradbasis(mesh);
femspace = get_femspace(mesh, vp.elemtype);
[vp, mesh] = init_poisson(mesh, femspace, vp);

%% ASSEMBLE STIFFNESS MATRIX (DIFFUSION)
A = assemble_stima(mesh, femspace, vp);
NA = size(A,1);

%% ASSEMBLE ADDITIONAL DG TERMS FOR DIFFUSION
%  1) FLUX TERMS
A_stima_dg_diff_flux = assemble_dg_diffusion_flux(mesh, femspace, vp);

if strcmp(vp.diffpart,'sipg')
    A = A - A_stima_dg_diff_flux - A_stima_dg_diff_flux.';
elseif strcmp(vp.diffpart,'nipg')
    A = A - A_stima_dg_diff_flux + A_stima_dg_diff_flux.';
elseif strcmp(vp.diffpart,'iipg')
    A = A - A_stima_dg_diff_flux;
else
    error('unknown diffusive dg part vp.diffpart');
end

%  2) PENALIZATION TERMS
A_stima_dg_diff_pen = assemble_dg_diffusion_pen(mesh, femspace, vp);
A = A + A_stima_dg_diff_pen;


%% ASSEMBLE RIGHT HAND SIDE
f = assemble_rhs(mesh, femspace, vp);

%% BOUNDARY CONDITIONS
poisson_boundary_conditions_dg;

%% LINEAR SYSTEM SOLUTION
ticsolver = tic;
switch vp.solver
	case 'backslash'
		sol = A \ permute(f, [1, 3, 2]);
		sol = permute(sol, [1, 3, 2]);
	case 'agmg'
		sol = zeros(size(A,1), 1, vp.NF);
		for i=1:size(f,3)
			sol(:,1,i) = agmg(A,f(:,1,i));
		end
end
if vp.ispureneumann % normalize
	sol = bsxfun(@minus,sol, get_average(mesh, femspace, sol));
end
stats.solver=toc(ticsolver);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function poisson_boundary_conditions_dg
        % FIND DIRICHLET DOF
        bdidx = zeros(NA,1);
        vp.ispureneumann = numel(mesh.bdflag{1}) == 0;
        if vp.ispureneumann
            bdidx(1) = 1;
        else
            [bddof, bdcoor] = get_bddof(mesh,femspace,mesh.bdflag{1});
            bdidx(bddof) = 1;
        end
        
        % (ALMOST EVERYWHERE) ZERO SOLUTION WITH RIGHT DIRICHLET BC
        uD = zeros(NA,1);
        if isnumeric(vp.bc{1})
            if vp.bc{1} ~= 0
                uD(bddof) = vp.bc{1};
            end
        else
            uD(bddof) = vp.bc{1}(bdcoor);
        end
        
        % UPDATE RHS BASED ON THE DIRICHLET BC
        f = bsxfun(@minus, f, A*uD);
        Tbd = spdiags(bdidx, 0, NA, NA);
        T = spdiags(1-bdidx, 0, NA, NA);
        A = T*A*T+ Tbd;
        
        % INSERT NEUMANN B.C.
        if ~isempty(mesh.bdflag{2})
            error('pure neumann B.C. not implemented for DG (not even for poisson)');
        end
        
        % ENFORCE STRONGLY THE DIRICHLET B.C.
        f(bdidx==1) = uD(bdidx==1);
    end

end