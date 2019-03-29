function res = get_form(mesh, whe, ...
	femspace1, sol1, der1, ...
	femspace2, sol2, der2)
%get_form gets a bilinear form
%% INITIALIZATION
dim = size(mesh.elem,2) - 1;
if ~isfield(mesh,'volume'), mesh.volume = simplex_volume(mesh); end
[lambda, weight] = quadpts(dim, femspace1.deg + femspace2.deg...
	-any(der1(:)~=0) -any(der2(:)~=0));

res = 0;
for i=1:numel(weight)
	res = res + weight(i) * ...
		evalf(mesh, whe, femspace1, lambda(i,:), sol1, der1) .* ...
		evalf(mesh, whe, femspace2, lambda(i,:), sol2, der2);
end
if isempty(whe)
	res = sum(bsxfun(@times, sum(res,2), mesh.volume), 1);
else
	res = sum(bsxfun(@times, sum(res,2), mesh.volume(whe)), 1);
end
end



