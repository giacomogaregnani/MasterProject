function fmingradp = get_fmingradp(mesh, vp, femspace, sol)
if ~isfield(mesh,'volume')
	mesh.volume = simplex_volume(mesh);
end
NQ  = numel(vp.aweight);
dim = size(mesh.elem,2) - 1;

fmingradp = 0;
for k=1:NQ
  ef = evalf(mesh, 'all', [], vp.alambda(k,:), vp.f);
  a = 0;
  for n=1:dim
    a = a + (ef(:,n+1) - evalf(mesh, 'all', femspace, ...
      vp.alambda(k,:), sol, n)).^2;
  end
  fmingradp = fmingradp + vp.aweight(k) * a;
end
fmingradp = fmingradp .* simplex_volume(mesh);
end

