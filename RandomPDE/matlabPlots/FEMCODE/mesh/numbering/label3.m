function mesh = label3(mesh)
%LABEL3 labels a 3D mesh in the sense of [AMP2000]
%
% mesh = label(mesh) renumbers the elements such that 
% mesh.elem(:,1) and mesh.elem(:,2) create the longest edge in each
% tetrahedron. It then adds three fields to the mesh structure:
% * flag: which is of size NT x 1 and set to false for every element
% * me1: which is of size NT x 1 and has value 1, 2, or 3, depending on
% which edge (with local indices) (1-3), (1-4), (3-4) is the longest.
% * me2: the same as me1 but with (2-3), (2,4), (3,4)
%
% [AMP2000] D. N. Arnold, A. Mukherjee, L. Pouly, Locally Adapted
% Tetrahedral Meshes Using Bisection, SIAM J. SCI. COMPUT., 2000, Vol. 22,
% No. 2, pp. 431--448

% ALREADY LABELED?
if isfield('mesh','flag'), return; end

%% PRE INITIALIZATION
% we assume that there is a lot of dummy elements (0,0,0,0) at the bottom
% of the array mesh.elem.
NTold = size(mesh.elem,1);
NT = find(~all(mesh.elem, 2), 1) - 1;
mesh.elem = mesh.elem(1:NT, :);

%% INITIALIZATION
dim = size(mesh.elem,2)-1;
%NT = size(mesh.elem,1);
pairs = dim*(dim+1)/2;

[edge, elem2edge] = auxstructure(mesh,'edge','elem2edge');
edgeLength = sum(get_vec(mesh,edge).^2,2);
[~, sortedEdges] = sort(edgeLength);
edgeLength(sortedEdges) = 1:size(edge,1);
[~, idx] = max(edgeLength(elem2edge),[],2);

regroup = [ get_subsimplices(1:dim+1, 1), ...
	flipdim( get_subsimplices(1:dim+1, dim - 2), 1)];
oddperm = permutationparity(regroup,2) == 1;
regroup(oddperm, [1 2]) =  regroup(oddperm, [2 1]);

for i=1:pairs % (can be run from 2, if orientation is correct for 1)
	mesh.elem(idx == i, 1:dim+1) = ...
		mesh.elem(idx == i, regroup(i,:));
end

% RECALCULATE LENGTHS
[edge, elem2edge] = auxstructure(mesh,'edge','elem2edge');
edgeLength = sum(get_vec(mesh,edge).^2,2);
[~, sortedEdges] = sort(edgeLength);
edgeLength(sortedEdges) = 1:size(edge,1);

[~, markedEdgeNonrefFace1] = max(edgeLength(elem2edge(:,[2,3,6])),[],2);
[~, markedEdgeNonrefFace2] = max(edgeLength(elem2edge(:,[4,5,6])),[],2);
mesh.flag = markedEdgeNonrefFace1 + 3*(markedEdgeNonrefFace2-1);
mesh.flag = uint8(mesh.flag);

% RESIZE BACK
if (NTold > NT)
	mesh.elem(NTold,dim+1) = 0;
end