function [Bq, Fq, micMesh, micVp, micFemspace] = get_forms(micMesh, micVp)

dim = size(micMesh.elem, 2) - 1;
micMesh = gradbasis(micMesh);
micFemspace = get_femspace(micMesh, micVp.elemtype);
[micVp, micMesh] = init_elasticity_orth(micMesh, micFemspace, micVp);
NA = micFemspace.ndof;
Bq = assemble_stiffness_2d(micMesh, micFemspace, micVp);
Fq = assemble_rhs_elasticity2d(micMesh, micFemspace, micVp);
Fq = reshape(Fq, [dim*NA, 1, micVp.NF]);

elasticity_boundary_conditions;

Fq = reshape(Fq, [dim*NA, micVp.NF]);

function elasticity_boundary_conditions
    micVp.ispureneumann = numel(micMesh.bdflag{1}) == 0;
    micVp.ispuredirichlet = numel(micMesh.bdflag{2}) == 0;
    if micVp.ispureneumann
        bddof = get_dof(micMesh, 1, micFemspace, 1);
        bdcoor = get_rc(micMesh, [], 1, [], micFemspace.nodelambda(1,:));
    else
        [bddof, bdcoor] = get_bddof(micMesh, micFemspace, micMesh.bdflag{1});
    end
    bdidx = zeros(dim*NA, 1);
    for n = 0:dim-1
        bdidx(NA*n + bddof) = 1;
    end
    uD = zeros(dim*NA, 1);
    % Insert Dirichlet values
    if ~micVp.ispureneumann
        if isnumeric(micVp.bc{1})
            for n = 1:dim
                uD((n-1)*NA + bddof) = micVp.bc{1}(n);
            end
        else
            dirval = micVp.bc{1}(bdcoor);
            for n = 1:dim
                uD((n-1)*NA + bddof) = dirval(:,n);
            end
        end
    end
    
    Fq = bsxfun(@minus, Fq, Bq*uD);
    Tbd = spdiags(bdidx, 0, dim*NA, dim*NA);
    T = spdiags(1-bdidx, 0, dim*NA, dim*NA);
    Bq = T*Bq*T + Tbd;
    if ~micVp.ispuredirichlet
        fi = zeros(dim*NA, 1);
        [lambda, weight] = quadpts(dim - 1, max(2*femspace.deg, 1));
        NQ = numel(weight);
        for j = 1:dim+1
            whe = micMesh.bdflag{2}(micMesh.bdflag{2}(:,2) == j, 1);
            if isempty(whe)
                continue;
            end
            lam = [lambda(:,1:j-1), zeros(size(lambda,1),1), lambda(:,j:dim)];
            ldof = find(micFemspace.nodelambda(:,j)==0)';
            ssvol = subsimplex_volume(micMesh, [], whe, [1:j-1, j+1:dim+1]);
            for k = 1:NQ
                if isnumeric(micVp.bc{2})
                    ef = micVp.bc{2};
                else
                    pxy = get_rc(micMesh, [], whe, [], lam(k,:));
                    if nargin(micVp.bc{2})==1
                        ef = micVp.bc{2}(pxy);
                    else
                        normals = get_normals(micMesh, whe, j);
                        ef = micVp.bc{2}(pxy, normals);
                    end
                end
                for n = 1:dim
                    if isnumeric(micVp.bc{2}) && (micVp.bc{2}(n)==0)
                        continue;
                    end
                    for m = ldof
                        gdof = get_dof(micMesh, whe, micFemspace, m);
                        fi((n-1)*NA + gdof) = fi((n-1)*NA + gdof) + ...
                            weight(k)*ssvol.*ef(:,n).*...
                            evalb(micMesh,whe,lam(k,:),m,0,micFemspace.elemtype);
                    end
                end
            end
            Fq = bsxfun(@plus, Fq, fi);
        end
    end
    % Enforce strongly Dirichlet
    for ii = 1:micVp.NF
        Fq(bdidx==1,:,ii) = uD(bdidx==1);
    end
end
end

