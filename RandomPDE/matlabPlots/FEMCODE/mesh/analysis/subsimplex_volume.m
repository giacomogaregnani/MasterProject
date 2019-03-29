function volume = subsimplex_volume(mesh, codim, whe, whi, sub)
%SUBSIMPLEX_VOLUME gives measures of subsimplices of a mesh
%
% volume = subsimplex_volume(mesh, codim) gives the measure of subsimplices
% of codimension codim in an NT x Nsub array. The subsimplices are ordered
% based on the lexicografical order of the local indices.
%
% volume = subsimplex_volume(mesh, 'opposite') takes codim = 1 but the
% order is reversed. 
%
% volume = subsimplex_volume(mesh, [], whe, whi) computes one volume for 
% each element in whe. The subsimplex is given by local node numbers: 
% a subset whi of 1:dim+1. whi can be 2-dimensional, specifying
% local node numbers for each element.
% 
% volume = subsimplex_volume(mesh, [], [], [], sub) computes one volume for
% each row of subs, which is a 2D array of global node indices.

if nargin < 2, codim = []; end
if nargin < 3, whe = []; end
if nargin < 4, whi = []; end
if nargin < 5, sub = []; end

d = size(mesh.node, 2);
NT = size(mesh.elem, 1);

if ~isempty(codim);
	if (codim == 0)
		volume = simplex_volume(mesh);
		return;
	end
	dim = size(mesh.elem, 2) - 1;
	if isnumeric(codim)
		dim2 = dim - codim;
		flip = false;
	else
		dim2 = dim - 1;
		flip = true;
	end	
	
	% ASSEMBLE VECTOR MATRIX
	s = get_subsimplices(1:dim+1, dim2);
	ns = size(s,1);
	if dim2 == 0
		volume = ones(NT,ns);
		return;
	end
	
	edge = get_subsimplices(1:dim+1, 1);
	nedge = size(edge, 1);
	
	vec = zeros(NT, d, nedge);
	for i=1:nedge
		vec(:,:,i) = get_vec(mesh, [], 'all', [edge(i,1), edge(i,2)]);
	end
	
	volume = zeros(NT, ns);
	for i=1:ns
		volume(:,i) = detn(vec(:,:,decode(s(i,1), s(i,2:end))));
	end
	
	volume = volume / factorial(dim2);
	if flip, 
		volume = flipdim(volume,2); 
	end
elseif ~isempty(sub) 
	NS = size(sub,1);
	dim2 = size(sub,2)-1;
	if dim2 == 0
		volume = ones(NS,1);
		return;
	end
	
	vec = zeros(NS, d, dim2);
	for i=1:dim2
		vec(:,:,i) = get_vec(mesh, sub(:,[1,i+1]));
	end
	volume = detn(vec);
	volume = volume / factorial(dim2);
else
	if isempty(whe) || strcmp(whe,'all')
		NS = NT;
	else
		NS = numel(whe);
	end
	dim2 = size(whi,2)-1;
	if dim2 == 0
		volume = ones(NS,1);
		return;
	end
	vec = zeros(NS, d, dim2);
	for i=1:dim2
		vec(:,:,i) = get_vec(mesh, [], whe, whi(:,[1,i+1]));
	end
	volume = detn(vec);
	volume = volume / factorial(dim2);
end		
		
%% DECODE FUNCTION
% dim+1 edges are sorted in edge: 1~2, 1~3, 1~4, ..., 1~dim+1, 2~3, 2~4,
% ..., 2~dim+1, 3~4, 3~5, ..., dim-1~dim, dim-1~dim+1, dim~dim+1. Suppose
% you have an edge a~b (or edges, in that case one or both are vectors).
% Function decode gives their index in the edge array.
function ind = decode(a,b)
ind = (dim+1)*(a-1) + (b-a) - a.*(a-1)/2;
end

end

