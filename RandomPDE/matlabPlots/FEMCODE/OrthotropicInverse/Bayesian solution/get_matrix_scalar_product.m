function W = get_matrix_scalar_product(mesh, vp)

dim = size(mesh.elem,2) - 1;
NT = size(mesh.elem,1);
unity_tensor = zeros(NT, dim+1, dim+1);
for k = 1:NT
    unity_tensor(k,:,:) = eye(dim+1);
end
femspace = get_femspace(mesh, vp.elemtype);
vp.aquad = max(2*(femspace.deg-1),1);
[vp.alambda, vp.aweight] = quadpts(dim, vp.aquad);
NQ = numel(vp.aweight);
vp.a = zeros(NT, dim+1, dim+1, NQ);
for qp = 1:NQ
    vp.a(:, :, :, qp) = unity_tensor;
end
W = get_forms(mesh, vp);