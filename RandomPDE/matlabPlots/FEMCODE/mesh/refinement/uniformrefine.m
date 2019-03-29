function [mesh, father] = uniformrefine(mesh, father)
%% UNIFORMREFINE  uniformly refine an n-D triangulation.
%
% [mesh] = uniformrefine(mesh) divides each triangle into 4 small
% similar triangles or each tetrahedron into eight smaller tetrahedrons
% (the division in not unique, we take the more regular).
%
% [mesh, bdFlag] = uniformrefine(mesh, bdFlag) updates boundary
% conditions represented by bdFlag.
%
% [mesh, ~, HB] = uniformrefine(mesh, [] , HB) outpus (and/or updates) HB 
% array which is useful for nodal interpolation and the father array.
%
% [mesh, ~, ~, father] = uniformrefine(mesh) outputs in addition the father 
% array.
%
% 1D refining scheme: 
%
% 1----3----2
%
% 2D refining scheme: 
%
%      3
%     / \ 
%    5---6
%   / \ / \
%  1---4---2
%
% 3D refining scheme: cut the four corner tetrahedra and then define the
% rest using the shortest diagonal
%
%          4--_
%         / \  10-__
%        /   \   __--3
%       7    _9--   /
%      /  _-6  \   8
%     /_--      \ /
%    1-----5-----2
%
% Reference: S. Zhang. Successive subdivisions of tetrahedra and multigrid
% methods on tetrahedral meshes. Houston J. Math. 21, 541?556, 1995.
%
% TODO:
% 1. see if one does not get a better and simpler algorithm using:
%
% Reference: J. Bey. Simplicial grid refinement: on Freudenthal's algorithm
% and the optimal number of congruence classes. Numer. Math.. 85(1):1--29,
% 2000. p11 Algorithm: RedRefinement3D.

mesh = cleanfields(mesh);
dim = size(mesh.elem,2)-1;
N = size(mesh.node,1); 
NT = size(mesh.elem,1); 

if nargin < 2, father = (1:NT)'; end
if nargout >= 2, father=repmat(father, [2^dim,1]); end

if dim > 1
	[sim1, elem2sim1] = auxstructure(mesh,'sim1','elem2sim1');
else
	sim1 = mesh.elem;
	elem2sim1 = (1:NT)';
end
NE = size(sim1,1);

%% Add new nodes
mesh.node(N+1:N+NE,:) = get_rc(mesh, sim1);
t = [mesh.elem, elem2sim1+N];

if dim == 1
	mesh.elem = [t(:,[1,3]); t(:,[3,2])];
elseif dim == 2
	mesh.elem = [t(:,[1,4,5]); t(:,[4,2,6]); t(:,[5,6,3]); t(:,[6,5,4])];
elseif dim == 3
	% 4 corner elements
	mesh.elem = ...
		[t(:,[1,5,6,7]); t(:,[5,2,8,9]); t(:,[6,8,3,10]); t(:,[7,9,10,4])];

	d1 = sum(get_vec(mesh, t(:,[7 8])).^2,2);
	d2 = sum(get_vec(mesh, t(:,[6 9])).^2,2);
	d3 = sum(get_vec(mesh, t(:,[5 10])).^2,2);
	[~,minindex] = min([d1 d2 d3],[],2);
		
	idx = find(minindex == 1); % diagonal 7-8 is the shortest
	mesh.elem(4*NT+idx,:) = t(idx,[7,8,10,9]); %1
	mesh.elem(5*NT+idx,:) = t(idx,[6,8,10,7]); %2
	mesh.elem(6*NT+idx,:) = t(idx,[9,7,8,5]); %3
	mesh.elem(7*NT+idx,:) = t(idx,[5,8,6,7]); %4
	
	idx = find(minindex == 2); % diagonal 6-9 is the shortest
	mesh.elem(4*NT+idx,:) = t(idx,[6,8,10,9]); %1
	mesh.elem(5*NT+idx,:) = t(idx,[6,9,10,7]); %2
	mesh.elem(6*NT+idx,:) = t(idx,[9,7,6,5]); %3
	mesh.elem(7*NT+idx,:) = t(idx,[5,8,6,9]); %4
	
	idx = find(minindex == 3); % diagonal 5-10 is the shortest
	mesh.elem(4*NT+idx,:) = t(idx,[5,8,10,9]); %1
	mesh.elem(5*NT+idx,:) = t(idx,[6,5,10,7]); %2
	mesh.elem(6*NT+idx,:) = t(idx,[9,7,10,5]); %3
	mesh.elem(7*NT+idx,:) = t(idx,[5,8,6,10]); %4
else
	error('DIMENSION > 3 NOT YET IMPLEMENTED');
end

% UPDATE BOUNDARY FLAG
if isfield(mesh,'bdflag') && ~ischar(mesh.bdflag)
	for j=1:numel(mesh.bdflag)
		whe = cell(dim+1,1);
		for k=1:dim+1
			whe{k} = mesh.bdflag{j}(mesh.bdflag{j}(:,2) == k, 1);
		end
		mesh.bdflag{j} = zeros(0,2);
		for k=1:dim+1
			for m=1:dim+1
				if k~=m
					mesh.bdflag{j} = [mesh.bdflag{j};
						whe{k} + (m-1)*NT, k*ones(size(whe{k},1),1)];
				end
			end
			if dim==3
				mesh.bdflag{j} = [mesh.bdflag{j};
					whe{k} + (4+k-1)*NT, k*ones(size(whe{k},1),1)];
			end
		end
	end
end
end
