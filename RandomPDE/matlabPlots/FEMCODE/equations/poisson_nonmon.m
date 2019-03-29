function [sol, femspace, mesh, vp] = poisson_nonmon(mesh, vp, options)
%POISSON solves a nonlinear nonmonotone problem given by mesh and vp

%% OPTIONS PARSING
if nargin<3, options = struct; end
if ~isfield(options,'verbose'), options.verbose = false; end
if ~isfield(options,'tolerance'), options.tolerance = 1e-6; end
if ~isfield(options,'maxIt'), options.maxIt = 100; end
dim = size(mesh.elem,2) - 1;

%% FINITE ELEMENTS INITIALIZATION
mesh = gradbasis(mesh);
femspace = get_femspace(mesh, vp.elemtype);
[vp, mesh] = init_poisson_nonmon(mesh, femspace, vp);

%% ASSEMBLE RIGHT HAND SIDE (INDEPENDENT OF NEWTON ITERATIONS)
rhs_f = assemble_rhs(mesh, femspace, vp);

%% NEWTON ITERATION
ND = femspace.ndof;
sol = zeros(ND,1);
if isfield(options,'sol') && size(options.sol,1) == ND, sol = options.sol(:,1); end

stats.increment = [];
increment = 1;
countIt = 0;


while increment > options.tolerance && countIt < options.maxIt
    
    countIt = countIt+1;
    solold = sol;
    
    %% ASSEMBLE STIFFNESS MATRIX
    A_stima = assemble_stima_nonmon(mesh, femspace, vp, solold);

    A_jac_stima = assemble_stima_jac_nonmon(mesh, femspace, vp, solold);
    
    A = A_stima + A_jac_stima;
    NA = size(A,1);
    f = A_stima*solold - rhs_f;
    
    %% BOUNDARY CONDITIONS
    elliptic_boundary_conditions;

    %% LINEAR SYSTEM SOLUTION
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
    
    increment = max(abs(sol));
    sol = solold - sol;
    stats.increment = [stats.increment, increment];
    
end

stats.newtonIt = countIt;


%% DISPLAY STATISTICS
if options.verbose
	fprintf('NEWTON IT:  %d, INCREMENT:  %1.4e, TOLERANCE:  %1.4e\n', ...
			stats.newtonIt, stats.increment(end), options.tolerance);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function elliptic_boundary_conditions
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
            error('neumann B.C. not implemented for elliptic problems');
        end
        
        % ENFORCE STRONGLY THE DIRICHLET B.C.
        f(bdidx==1) = uD(bdidx==1);
    end

end