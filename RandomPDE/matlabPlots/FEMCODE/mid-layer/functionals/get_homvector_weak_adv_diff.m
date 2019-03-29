function b0vector = get_homvector_weak_adv_diff(sol, femspace, mesh, vp)
% Calculates homogenized tensor

dim = size(mesh.elem,2) - 1;
b0vector = zeros(dim,1);

%% CALCULATE HOMOGENIZED TENSOR
area = sum(mesh.volume);
for i=1:dim
    res = 0;
    for s=1:numel(vp.aweight)
        res = res + vp.aweight(s) * ...
            evalf(mesh, 'all', femspace, vp.alambda(s,:), vp.b, [], i);
        for m = 1:dim
            res = res + vp.aweight(s) * ...
                evalf(mesh, 'all', femspace, vp.alambda(s,:), vp.b, [], m) .* ...
                evalf(mesh, 'all', femspace, vp.alambda(s,:), sol(:,:,i), m);
        end
    end
    b0vector(i) = dot(res, mesh.volume)/area;
end

end