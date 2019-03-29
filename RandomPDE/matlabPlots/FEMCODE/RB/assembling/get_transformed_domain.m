function mesh = get_transformed_domain(RB, par)
%GET_TRANSFORMED_DOMAIN Summary of this function goes here
%   Detailed explanation goes here

%% better version
[mesh, father] = RB.getfine(RB);
R = numel(RB.C);
dim = size(mesh.elem,2)-1;
done = [];
for i=1:R
  C = subsmu(RB.C{i}, par);
  G = subsmu(RB.G{i}, par);
  for j=1:dim+1
    el = find(father==i);
    nodes = mesh.elem(el, j);
    [newnodes, ind] = setdiff(nodes,done);
    done = union(done, newnodes);
    nodeval = get_rc(mesh, [], el(ind), j, 1, struct('adjust_to','bary'));
    mesh.node(newnodes,:) = bsxfun(@plus, C', nodeval * G');
  end
end

end