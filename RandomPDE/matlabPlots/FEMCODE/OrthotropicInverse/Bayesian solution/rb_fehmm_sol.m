function macSol = rb_fehmm_sol(macMesh, macVp, parameter, AqRB, FqRB, Q, micMesh, micVp, micFemspace, RBspace)

dim = size(macMesh.elem, 2) - 1;
macFemspace = get_femspace(macMesh, macVp.elemtype);
macVp.aquad = max(2*(macFemspace.deg-1),1);
[macVp.alambda, macVp.aweight] = quadpts(dim, macVp.aquad);
macVp.fully_discrete = true;

NQ = numel(macVp.aweight);
NT = size(macMesh.elem, 1);
NR = size(AqRB,1);
NdofMic = micFemspace.ndof;
res = zeros(NT, Q);
for q = 1:Q
    for qp = 1:NQ
        res(:,q) = res(:,q) + macVp.aweight(qp) * (evalf(macMesh, 'all', macFemspace, macVp.alambda(qp,:), micVp.ThetaField(macMesh.node, q, parameter)));
    end
end
% Change no dependence on space
xloc = get_rc(macMesh, [], 'all', [], macVp.alambda(qp,:));
a0tensor = zeros(NT, dim+1, dim+1);
ARB = reshape(sum(bsxfun(@times, res(1,:)', permute(AqRB, [3, 1, 2])), 1), NR, NR);
FRB = reshape(sum(bsxfun(@times, res(1,:)', permute(FqRB, [3, 1, 2])), 1), NR, 1, micVp.NF);
micSolX = ARB\permute(FRB, [1 3 2]);
micSolX = RBspace*micSolX;
micSolX = reshape(micSolX, NdofMic, dim, []);
micVpX = micVp;
micVpX.a = @(x, kk, ll) feval(micVp.a_full, xloc(1,:), x, kk, ll, parameter);
a0tensor(1, :, :) = get_homogenized_tensor_elasticity(micSolX, micFemspace, micMesh, micVpX);
for k = 2:NT
    a0tensor(k,:,:) = a0tensor(1,:,:);
end
macVp.a(:,:,:,1) = a0tensor;

% for qp = 1:NQ
%     xloc = get_rc(macMesh, [], 'all', [], macVp.alambda(qp,:));
%     a0tensor = zeros(NT, dim+1, dim+1);
%     for k = 1:NT
%         fprintf('QP: %.1d / %.1d;  Elem: %.5d / %.5d\n', ...
%             qp, NQ, k, NT);
%         ARB = reshape(sum(bsxfun(@times, res(k,:)', permute(AqRB, [3, 1, 2])), 1), NR, NR);
%         FRB = reshape(sum(bsxfun(@times, res(k,:)', permute(FqRB, [3, 1, 2])), 1), NR, 1, micVp.NF);
%         micSolX = ARB\permute(FRB, [1 3 2]);
%         micSolX = RBspace*micSolX;
%         micSolX = reshape(micSolX, NdofMic, dim, []);
%         micVpX = micVp;
%         micVpX.a = @(x, kk, ll) feval(micVp.a_full, xloc(k,:), x, kk, ll, parameter);
%         a0tensor(k, :, :) = get_homogenized_tensor_elasticity(micSolX, micFemspace, micMesh, micVpX);
%     end
%     macVp.a(:,:,:,qp) = a0tensor;
% end

[macSol, macFemspace, macMesh, macVp] = my_elasticity_2d(macMesh, macVp);

end