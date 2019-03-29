function [macSol, macFemspace, macMesh, macVp, a0tensor] = fehmm_elasticity(macMesh, macVp, micMesh, micVp)

dim = size(macMesh.elem, 2) - 1;
macFemspace = get_femspace(macMesh, macVp.elemtype);
macVp.aquad = max(2*(macFemspace.deg-1),1);
[macVp.alambda, macVp.aweight] = quadpts(dim, macVp.aquad);
macVp.fully_discrete = true;

NQ = numel(macVp.aweight);
NT = size(macMesh.elem, 1);

for qp= 1:NQ
    
    xloc = get_rc(macMesh, [], 'all', [], macVp.alambda(qp,:));
    a0tensor = zeros(NT, dim + 1, dim + 1); % in two dimensions the 4th order tensor can be rewritten as 3x3 matrix (9 entries)
    
    for k = 1:NT
        fprintf('QP: %.1d / %.1d;  Elem: %.5d / %.5d\n', ...
     qp, NQ, k, NT);
 
        micVpX = micVp;
        micVpX.a = @(x, kk, ll) feval(micVp.a, xloc(k,:), x, kk, ll);
        [micSolX, micFemspaceX, micMeshX, micVpX] = micro_problem_elasticity(micMesh, micVpX); % --> write this function
        a0tensor(k,:,:) = get_homogenized_tensor_elasticity(micSolX, micFemspaceX, micMeshX, micVpX); % --> write this function
    end

    macVp.a(:, :, :, qp) = a0tensor;
    
end

[macSol, macFemspace, macMesh, macVp] = my_elasticity_2d(macMesh, macVp); % --> allow fully discrete tensor

end