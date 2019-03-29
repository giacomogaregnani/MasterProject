function a0tensor = get_homogenized_tensor_elasticity(sol, femspace, mesh, vp)

dim = size(mesh.elem,2) - 1;
a0tensor = zeros(dim+1, dim+1);
area = sum(mesh.volume);
area = 1;
% for i = 1:dim
%     for j = 1:dim
%         res = 0;
%         for s = 1:numel(vp.aweight)
%             res = res + vp.aweight(s) * evalf(mesh, 'all', femspace, vp.alambda(s,:), vp.a, [], i, j);
%             for m = 1:dim
%                 res = res + vp.aweight(s) * evalf(mesh, 'all', femspace, vp.alambda(s,:), vp.a, [], i, m) .* ...
%                     evalf(mesh, 'all', femspace, vp.alambda(s,:), sol(:,m,j), m);
%             end
%         end
%         a0tensor(i,j) = dot(res, mesh.volume)/area;
%     end
% end
res = 0;
for s = 1:numel(vp.aweight)
    res = res + vp.aweight(s) * evalf(mesh, 'all', femspace, vp.alambda(s,:), vp.a, [], 1, 1);
    res = res + vp.aweight(s) * evalf(mesh, 'all', femspace, vp.alambda(s,:), vp.a, [], 1, 1) .* ...
        evalf(mesh, 'all', femspace, vp.alambda(s,:), sol(:,1,1), 1);
    res = res + vp.aweight(s) * evalf(mesh, 'all', femspace, vp.alambda(s,:), vp.a, [], 1, 2) .* ...
        evalf(mesh, 'all', femspace, vp.alambda(s,:), sol(:,2,1), 2);
end
a0tensor(1,1) = dot(res, mesh.volume)/area;

res = 0;
for s = 1:numel(vp.aweight)
    res = res + vp.aweight(s) * evalf(mesh, 'all', femspace, vp.alambda(s,:), vp.a, [], 1, 2);
    res = res + vp.aweight(s) * evalf(mesh, 'all', femspace, vp.alambda(s,:), vp.a, [], 1, 1) .* ...
        evalf(mesh, 'all', femspace, vp.alambda(s,:), sol(:,1,2), 1);
    res = res + vp.aweight(s) * evalf(mesh, 'all', femspace, vp.alambda(s,:), vp.a, [], 1, 2) .* ...
        evalf(mesh, 'all', femspace, vp.alambda(s,:), sol(:,2,2), 2);
end
a0tensor(1,2) = dot(res, mesh.volume)/area;

res = 0;
for s = 1:numel(vp.aweight)
    res = res + vp.aweight(s) * evalf(mesh, 'all', femspace, vp.alambda(s,:), vp.a, [], 2, 1);
    res = res + vp.aweight(s) * evalf(mesh, 'all', femspace, vp.alambda(s,:), vp.a, [], 2, 1) .* ...
        evalf(mesh, 'all', femspace, vp.alambda(s,:), sol(:,1,1), 1);
    res = res + vp.aweight(s) * evalf(mesh, 'all', femspace, vp.alambda(s,:), vp.a, [], 2, 2) .* ...
        evalf(mesh, 'all', femspace, vp.alambda(s,:), sol(:,2,1), 2);
end
a0tensor(2,1) = dot(res, mesh.volume)/area;

res = 0;
for s = 1:numel(vp.aweight)
    res = res + vp.aweight(s) * evalf(mesh, 'all', femspace, vp.alambda(s,:), vp.a, [], 2, 2);
    res = res + vp.aweight(s) * evalf(mesh, 'all', femspace, vp.alambda(s,:), vp.a, [], 2, 1) .* ...
        evalf(mesh, 'all', femspace, vp.alambda(s,:), sol(:,1,2), 1);
    res = res + vp.aweight(s) * evalf(mesh, 'all', femspace, vp.alambda(s,:), vp.a, [], 2, 2) .* ...
        evalf(mesh, 'all', femspace, vp.alambda(s,:), sol(:,2,2), 2);
end
a0tensor(2,2) = dot(res, mesh.volume)/area;

res = 0;
for s = 1:numel(vp.aweight)
    res = res + vp.aweight(s) * evalf(mesh, 'all', femspace, vp.alambda(s,:), vp.a, [], 3, 3);
    res = res + vp.aweight(s) * evalf(mesh, 'all', femspace, vp.alambda(s,:), vp.a, [], 3, 3) .* ...
        evalf(mesh, 'all', femspace, vp.alambda(s,:), sol(:,1,3), 2);
    res = res + vp.aweight(s) * evalf(mesh, 'all', femspace, vp.alambda(s,:), vp.a, [], 3, 3) .* ...
        evalf(mesh, 'all', femspace, vp.alambda(s,:), sol(:,2,3), 1);
end
a0tensor(3,3) = dot(res, mesh.volume)/area;

end