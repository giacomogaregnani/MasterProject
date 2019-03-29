function vec = get_vec(mesh, sub, whe, whi)
%GET_VEC returns vectors defined by pairs of nodes 
%
% vec = get_vec(mesh, sub) computes the vectors starting at nodes sub(:,1) 
% and ending in sub(:,2)
%
% vec = get_vec(mesh, [], whe, whi) computes the vectors in elements whe
% between nodes with local numbers whi( ,1) and whi( ,2)
%
% CONVENTION: whe = [] implies whe = 1:NT
%
%   if numel(whi) == 2, then vec(i,:) is the vector joining 
%   mesh.elem(whe(i), whi(1)) and mesh.elem(whe(i), whi(2)). Otherwise, it 
%   is the vector joining mesh.elem(whe(i), whi(i, 1)) and 
%   mesh.elem(whe(i), whi(i, 2)).

if (nargin >= 2) && ~isempty(sub)
	vec = mesh.node(sub(:,2),:) - mesh.node(sub(:,1),:);
else
	NT= size(mesh.elem,1);
	if numel(whi) == 2
		if ~strcmp(whe,'all')
			vec = mesh.node(mesh.elem(whe,whi(2)),:) - ...
				mesh.node(mesh.elem(whe,whi(1)),:);
		else
			vec = mesh.node(mesh.elem(:,whi(2)),:) - ...
				mesh.node(mesh.elem(:,whi(1)),:);
		end
	else
		if ~strcmp(whe,'all')
			vec = mesh.node(mesh.elem( whe+NT*(whi(:,2)-1) ),:) - ...
				mesh.node(mesh.elem( whe+NT*(whi(:,1)-1) ),:);
		else
			vec = mesh.node(mesh.elem( (1:NT)'+NT*(whi(:,2)-1) ),:) - ...
				mesh.node(mesh.elem( (1:NT)'+NT*(whi(:,1)-1) ),:);
		end
	end
end
if isfield(mesh,'periodic') && mesh.periodic
	vec = per_diff(vec, mesh.box);
end
end

