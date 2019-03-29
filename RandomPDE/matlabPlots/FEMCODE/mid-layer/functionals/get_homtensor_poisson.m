function a0tensor = get_homtensor_poisson(sol, femspace, mesh, vp)
% Calculates homogenized tensor

dim = size(mesh.elem,2) - 1;
a0tensor = zeros(dim, dim);

%% CALCULATE HOMOGENIZED TENSOR
area = sum(mesh.volume);
for i=1:dim
    for j=1:dim
        res = 0;
        for s=1:numel(vp.aweight)
            res = res + vp.aweight(s) * ...
                evalf(mesh, 'all', femspace, vp.alambda(s,:), vp.a, [], i, j);
            for m = 1:dim
                res = res + vp.aweight(s) * ...
                    evalf(mesh, 'all', femspace, vp.alambda(s,:), vp.a, [], i, m) .* ...
                    evalf(mesh, 'all', femspace, vp.alambda(s,:), sol(:,:,j), m);
            end
        end
        a0tensor(i,j) = dot(res, mesh.volume)/area;
    end
end

end