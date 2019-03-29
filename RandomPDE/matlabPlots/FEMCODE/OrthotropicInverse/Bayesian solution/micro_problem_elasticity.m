function [sol, femspace, mesh, vp] = micro_problem_elasticity(mesh, vp)

dim = size(mesh.elem, 2) - 1;
mesh = gradbasis(mesh);
femspace = get_femspace(mesh, vp.elemtype);
[vp, mesh] = init_elasticity_orth(mesh, femspace, vp);
NA = femspace.ndof;
A = assemble_stiffness_2d(mesh, femspace, vp);
f = assemble_rhs_elasticity2d(mesh, femspace, vp);
f = reshape(f, [dim*NA, 1, vp.NF]);

elasticity_boundary_conditions;

sol = A\permute(f, [1 3 2]);
sol = reshape(sol, NA, dim, []);

if vp.ispureneumann
  sol = bsxfun(@minus,sol, get_average(mesh, femspace, sol));
end

function elasticity_boundary_conditions
    vp.ispureneumann = numel(mesh.bdflag{1}) == 0;
    vp.ispuredirichlet = numel(mesh.bdflag{2}) == 0;
    if vp.ispureneumann
        bddof = get_dof(mesh, 1, femspace, 1);
        bdcoor = get_rc(mesh, [], 1, [], femspace.nodelambda(1,:));
    else
        [bddof, bdcoor] = get_bddof(mesh, femspace, mesh.bdflag{1});
    end
    bdidx = zeros(dim*NA, 1);
    for n = 0:dim-1
        bdidx(NA*n + bddof) = 1;
    end
    uD = zeros(dim*NA, 1);
    % Insert Dirichlet values
    if ~vp.ispureneumann
        if isnumeric(vp.bc{1})
            for n = 1:dim
                uD((n-1)*NA + bddof) = vp.bc{1}(n);
            end
        else
            dirval = vp.bc{1}(bdcoor);
            for n = 1:dim
                uD((n-1)*NA + bddof) = dirval(:,n);
            end
        end
    end
    
    f = bsxfun(@minus, f, A*uD);
    Tbd = spdiags(bdidx, 0, dim*NA, dim*NA);
    T = spdiags(1-bdidx, 0, dim*NA, dim*NA);
    A = T*A*T + Tbd;
    if ~vp.ispuredirichlet
        fi = zeros(dim*NA, 1);
        [lambda, weight] = quadpts(dim - 1, max(2*femspace.deg, 1));
        NQ = numel(weight);
        for j = 1:dim+1
            whe = mesh.bdflag{2}(mesh.bdflag{2}(:,2) == j, 1);
            if isempty(whe)
                continue;
            end
            lam = [lambda(:,1:j-1), zeros(size(lambda,1),1), lambda(:,j:dim)];
            ldof = find(femspace.nodelambda(:,j)==0)';
            ssvol = subsimplex_volume(mesh, [], whe, [1:j-1, j+1:dim+1]);
            for k = 1:NQ
                if isnumeric(vp.bc{2})
                    ef = vp.bc{2};
                else
                    pxy = get_rc(mesh, [], whe, [], lam(k,:));
                    if nargin(vp.bc{2})==1
                        ef = vp.bc{2}(pxy);
                    else
                        normals = get_normals(mesh, whe, j);
                        ef = vp.bc{2}(pxy, normals);
                    end
                end
                for n = 1:dim
                    if isnumeric(vp.bc{2}) && (vp.bc{2}(n)==0)
                        continue;
                    end
                    for m = ldof
                        gdof = get_dof(mesh, whe, femspace, m);
                        fi((n-1)*NA + gdof) = fi((n-1)*NA + gdof) + ...
                            weight(k)*ssvol.*ef(:,n).*...
                            evalb(mesh,whe,lam(k,:),m,0,femspace.elemtype);
                    end
                end
            end
            f = bsxfun(@plus, f, fi);
        end
    end
    % Enforce strongly Dirichlet
    for ii = 1:vp.NF
        f(bdidx==1,:,ii) = uD(bdidx==1);
    end
end
end