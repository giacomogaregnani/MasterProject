function [bddof, bdcoor] = get_bddof(mesh, femspace, flag)
%GET_BDDOF Summary of this function goes here
%   Detailed explanation goes here

dim = size(mesh.elem,2) - 1;
d = size(mesh.node,2);
bddof = zeros(0,1);
bdcoor = zeros(0,d);
for i=1:dim+1
	f = flag(flag(:,2) == i, 1);
	if isempty(f)
		continue;
	end
	ldof = find(femspace.nodelambda(:,i)==0)';
	for j=ldof;
		bddof = [bddof; get_dof(mesh, f(:,1), femspace, j)];
		bdcoor = [bdcoor; get_rc(mesh, [], f(:,1), [], ...
			femspace.nodelambda(j,:))];
	end
end
[bddof,ind] = unique(bddof,'legacy');
bdcoor = bdcoor(ind,:);
end

