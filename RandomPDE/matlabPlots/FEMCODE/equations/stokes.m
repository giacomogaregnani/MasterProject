function [usol, ufemspace, psol, pfemspace, mesh, vp] = ...
	stokes(mesh, vp, options)
%STOKES solves a Stokes problem given by mesh and vp
%   Detailed explanation goes here

%% OPTIONS PARSING
if nargin<3, options = struct; end
if ~isfield(options,'verbose'), options.verbose = false; end

%% FINITE ELEMENTS INITIALIZATION
mesh = gradbasis(mesh);
ufemspace = get_femspace(mesh, vp.elemtype);
pfemspace = get_femspace(mesh, vp.pelemtype);
[vp, mesh] = init_stokes(mesh, ufemspace, vp);

%% CONSTANTS
dim = size(mesh.elem, 2) -1;
NA = ufemspace.ndof;
NB = pfemspace.ndof; 

%% ASSEMBLE STIFFNESS MATRIX
A = assemble_stima(mesh, ufemspace, vp);
A = stack_blkdiag(A,dim);

%% ASSEMBLE RIGHT HAND SIDE
f = assemble_rhs(mesh, ufemspace, vp);
f = reshape(f, [dim*NA, 1, vp.NF]);

%% ASSEMBLE DIV MATRIX
B = assemble_div(mesh, ufemspace, pfemspace, vp); 
g = zeros(NB, 1, vp.NF);

%% BOUNDARY CONDITIONS
stokes_boundary_conditions;

%% LINEAR SYSTEM SOLUTION
%if size(A,1) > 400000 && all(g(:)==0)
%	vp.solver = 'iuzawapcg';
%end
if ~all(g(:)==0)
	vp.solver = 'backslash';
end

switch vp.solver
	case 'backslash'
		A = [A, B'; B, sparse(NB, NB)];
		rhs = cat(1,f,g);
		if vp.ispuredirichlet
			solution = zeros(size(A,1), vp.NF);
			solution(1:end-1,:) = A(1:end-1,1:end-1) \ rhs(1:end-1,:);
		else
			solution = A \ rhs;
		end
		usol = reshape(solution(1:dim*NA,:), NA, dim, []);
		psol = reshape(solution(dim*NA+1:end,:), NB, 1, []);
		psol = bsxfun(@minus, psol, get_average(mesh, pfemspace, psol));
	case 'iuzawapcg'
		M = assemble_mass(mesh, pfemspace);
		M = spdiags(diag(M),0,size(B,1), size(B,1));
		if ~isfield(vp,'iuzawapcgopt'), 
			vp.iuzawapcgopt = struct; 
		end
		[usol, psol] = iuzawapcg(A, B, f, M, vp.iuzawapcgopt);
		usol = reshape(usol, NA, dim, []);
		psol = reshape(psol, NB, 1, []);
		psol = bsxfun(@minus, psol, get_average(mesh, pfemspace, psol));
	otherwise
		error('Unknown linear algebra solver');
end

%% NORMALIZE VELOCITY IN PURE NUEMANN PROBLEM
if any(vp.ispureneumann)
	usol(:,vp.ispureneumann,:) = ...
        bsxfun(@minus, usol(:,vp.ispureneumann,:), ...
		get_average(mesh, ufemspace, usol(:,vp.ispureneumann,:)));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	function stokes_boundary_conditions
	%% BOUNDARY CONDITIONS
	% FIND DIRICHLET DOF
    vp.ispureneumann = numel(mesh.bdflag{1}) == 0;
    vp.ispuredirichlet = numel(mesh.bdflag{2}) == 0;
    if vp.ispureneumann
        bddof = get_dof(mesh, 1, ufemspace, 1);
        bdcoor = get_rc(mesh, [], 1, [], ufemspace.nodelambda(1,:));
    else
        [bddof, bdcoor] = get_bddof(mesh, ufemspace, mesh.bdflag{1});
    end
    
    % boundary "logical" field
    bdidx = zeros(dim*NA,1);
    for n=0:dim-1
        bdidx(NA*n + bddof) = 1;
    end
    
    % insert Dirichlet values
    uD = zeros(dim*NA,1);
    if ~vp.ispureneumann
        if isnumeric(vp.bc{1})
            for n=1:dim
                uD((n-1)*NA + bddof) = vp.bc{1}(n);
            end
        else
            dirval = vp.bc{1}(bdcoor);
            for n=1:dim
                uD((n-1)*NA + bddof) = dirval(:,n);
            end
        end
    end
	
	% UPDATE RHS BASED ON THE DIRICHLET BC
	f = bsxfun(@minus,f,A*uD);
	g = bsxfun(@minus,g,B*uD);
	
	Tbd = spdiags(bdidx, 0, dim*NA, dim*NA);
	T = spdiags(1-bdidx, 0, dim*NA, dim*NA);
	A = T*A*T+ Tbd;
	B = B*T;
	
    % INSERT NEUMANN B.C.
    if ~vp.ispuredirichlet
        fi = zeros(dim*NA,1);
        [lambda, weight] = quadpts(dim-1, max(2*ufemspace.deg,1));
        NQ = numel(weight);
        for j=1:dim+1
            whe = mesh.bdflag{2}(mesh.bdflag{2}(:,2) == j, 1);
            if isempty(whe),
                continue;
            end
            lam = [lambda(:,1:j-1), zeros(size(lambda,1), 1), ...
                lambda(:,j:dim)];
            ldof = find(ufemspace.nodelambda(:,j)==0)';
            ssvol = subsimplex_volume(mesh, [], whe, [1:j-1, j+1:dim+1]);            
            for k=1:NQ
                if isnumeric(vp.bc{2})
                    ef = vp.bc{2};
                else
                    pxy = get_rc(mesh, [], whe, [], lam(k,:));
                    if nargin(vp.bc{2}) == 1
                        ef = vp.bc{2}(pxy);
                    else
                        normals = get_normals(mesh, whe, j);
                        ef = vp.bc{2}(pxy,normals);
                    end
                end
                for n=1:dim
                    if isnumeric(vp.bc{2}) && (vp.bc{2}(n) == 0),
                        continue;
                    end   
                    for m=ldof
                        gdof = get_dof(mesh, whe, ufemspace, m);
                        fi((n-1)*NA + gdof) = fi((n-1)*NA + gdof) + ...
                            weight(k) * ssvol .* ef(:,n) .*...
                            evalb(mesh, whe, lam(k,:), m, 0, ...
                            ufemspace.elemtype);
                    end
                end
            end            
        end
        f = bsxfun(@plus, f, fi);
    end
    
    % ENFORCE STRONGLY DIRICHLET B.C.
    for ii=1:vp.NF
      f(bdidx==1,:,ii) = uD(bdidx==1);
    end
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

