function [sol, femspace, mesh, vp, A, f] = poisson(mesh, vp, options)
%POISSON solves a Poisson problem given by mesh and vp

%% OPTIONS PARSING
if nargin<3, options = struct; end
if ~isfield(options,'verbose'), options.verbose = false; end
dim = size(mesh.elem,2) - 1;

%% FINITE ELEMENTS INITIALIZATION
mesh = gradbasis(mesh);
femspace = get_femspace(mesh, vp.elemtype);
[vp, mesh] = init_poisson(mesh, femspace, vp);

%% ASSEMBLE STIFFNESS MATRIX AND RHS
A = assemble_stima(mesh, femspace, vp);
f = assemble_rhs(mesh, femspace, vp);
NA = size(A,1);

%% BOUNDARY CONDITIONS
poisson_boundary_conditions;

%% LINEAR SYSTEM SOLUTION
switch vp.solver
  case 'backslash'
    sol = A \ permute(f, [1, 3, 2]);
    sol = permute(sol, [1, 3, 2]);
  case 'agmg'
    sol = zeros(NA, 1, vp.NF);
    for i=1:size(f,3)
      sol(:,1,i) = agmg(A,f(:,1,i));
    end
end

%% NORMALIZE SOLUTION IN PURE NUEMANN PROBLEM
if vp.ispureneumann % normalize
  sol = bsxfun(@minus,sol, get_average(mesh, femspace, sol));
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  function poisson_boundary_conditions
    % FIND DIRICHLET DOF
    vp.ispureneumann = numel(mesh.bdflag{1}) == 0;
    if vp.ispureneumann
      bddof = get_dof(mesh, 1, femspace, 1);
      bdcoor = get_rc(mesh, [], 1, [], femspace.nodelambda(1,:));
    else
      [bddof, bdcoor] = get_bddof(mesh,femspace,mesh.bdflag{1});
    end
    
    % boundary "logical" field
    bdidx = zeros(NA,1);
    bdidx(bddof) = 1;
    
    % (ALMOST EVERYWHERE) ZERO SOLUTION WITH RIGHT DIRICHLET BC
    uD = zeros(NA,1);
    if ~vp.ispureneumann
      if isnumeric(vp.bc{1})
        if vp.bc{1} ~= 0
          uD(bddof) = vp.bc{1};
        end
      else
        uD(bddof) = vp.bc{1}(bdcoor);
      end
    end
    
    % UPDATE RHS BASED ON THE DIRICHLET BC
    f = bsxfun(@minus, f, A*uD);
    Tbd = spdiags(bdidx, 0, NA, NA);
    T = spdiags(1-bdidx, 0, NA, NA);
    A = T*A*T+ Tbd;
    
    % INSERT NEUMANN B.C.
    [lambda, weight] = quadpts(dim-1, max(2*femspace.deg,1));
    NQ = numel(weight);
    if ~(isnumeric(vp.bc{2}) && (vp.bc{2} == 0))
      fi = zeros(NA,1);
      for j=1:dim+1
        whe = mesh.bdflag{2}(mesh.bdflag{2}(:,2) == j, 1);
        if isempty(whe),
          continue;
        end
        lam = [lambda(:,1:j-1), zeros(size(lambda,1),1), lambda(:,j:dim)];
        ldof = find(femspace.nodelambda(:,j)==0)';
        ssvol = subsimplex_volume(mesh, [], whe, [1:j-1,j+1:dim+1]);
        for k=1:NQ
          if isnumeric(vp.bc{2})
            ef = repmat(vp.bc{2}, size(whe));
          else
            pxy = get_rc(mesh, [], whe, [], lam(k,:));
            if nargin(vp.bc{2}) == 1
              ef = vp.bc{2}(pxy);
            else
              normals = get_normals(mesh, whe, j);
              ef = vp.bc{2}(pxy,normals);
            end
          end
          for m=ldof
            dof = get_dof(mesh, whe, femspace, m);
            fi(dof) = fi(dof) + weight(k) * ssvol .* ef .*...
              evalb(mesh, whe, lam(k,:), m, 0, femspace.elemtype);
          end
        end
      end
      f = bsxfun(@plus, f, fi);
    end
    
    % ENFORCE STRONGLY THE DIRICHLET B.C.
    for ii=1:vp.NF
      f(bdidx==1,:,ii) = uD(bdidx==1);
    end
  end
end