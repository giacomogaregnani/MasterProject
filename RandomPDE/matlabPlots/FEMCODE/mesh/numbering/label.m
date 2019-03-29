function [mesh, bdFlag] = label(mesh, whe, bdFlag)
%LABEL relabels elements such that first two nodes create the longest edge
%
% mesh = label(mesh) renumbers the elements such that 
% mesh.elem(:,1) and mesh.elem(:,2) create the longest edge in each simplex
%
% [mesh] = label(mesh,whe) does the same but only for a subset of elements
% specified in whe.
%  
% [mesh, bdFlag] = label(mesh, whe, bdFlag) also updates the array bdFlag


%% INITIALIZATION
dim = size(mesh.elem,2)-1;
NT = size(mesh.elem,1);
pairs = dim*(dim+1)/2;

if (dim == 1), return; end
if (nargin < 2) || strcmp(whe,'all'), whe = 1:NT; end
if islogical(whe), whe = find(whe); end
if (nargin<3), bdFlag = []; end

[edge, elem2edge] = auxstructure(mesh,'edge','elem2edge');
edgeLength = sum(get_vec(mesh,edge).^2,2);
[~, elemEdgeLength] = sort(edgeLength);
elemEdgeLength(elemEdgeLength) = 1:size(edge,1);
elemEdgeLength = elemEdgeLength(elem2edge);

% allEdge = get_subsimplices(mesh.elem(whe,:), 1);
% elemEdgeLength = reshape(sum(get_vec(mesh,allEdge).^2,2), NT, pairs);
[~, idx] = max(elemEdgeLength,[],2);

regroup = [ get_subsimplices(1:dim+1, 1), ...
	flipdim( get_subsimplices(1:dim+1, dim - 2), 1)];
oddperm = permutationparity(regroup,2) == 1;
regroup(oddperm, [1 2]) =  regroup(oddperm, [2 1]);

for i=1:pairs % (can be run from 2, if orientation is correct for 1)
	mesh.elem(whe(idx == i), 1:dim+1) = ...
		mesh.elem(whe(idx == i), regroup(i,:));
end
if ~isempty(bdFlag)
	for i=1:size(pairs,1)
		bdFlag(whe(idx == i), 1:dim+1) = ...
			bdFlag(whe(idx == i), regroup(i,:));
	end
end

% if (dim == 3)
%     % RECALCULATE LENGTHS
%     [edge, elem2edge] = auxstructure(mesh,'edge','elem2edge');
%     edgeLength = sum(get_vec(mesh,edge).^2,2);
%     [~,elemEdgeLength] = sort(edgeLength);
%     elemEdgeLength = elemEdgeLength(elem2edge);
%     
%     [~, mesh.me1] = max(elemEdgeLength(:,[2,3,6]),[],2);
%     [~, mesh.me2] = max(elemEdgeLength(:,[4,5,6]),[],2);
%     mesh.flag = false(NT,1);    
% end

