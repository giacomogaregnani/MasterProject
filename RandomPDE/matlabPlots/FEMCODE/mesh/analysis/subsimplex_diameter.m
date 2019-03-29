function diam = subsimplex_diameter(mesh, codim, ~, ~, sub)
%SUBSIMPLEX_VOLUME gives measures of subsimplices of a mesh
%
% volume = subsimplex_diameter(mesh, codim) gives the measure of subsimplices
% of codimension codim in an NT x Nsub array. The subsimplices are ordered
% based on the lexicografical order of the local indices.
%
% volume = subsimplex_diameter(mesh, 'opposite') takes codim = 1 but the
% order is reversed. 


if nargin < 2, codim = 'opposite'; end
if nargin < 5, sub = []; end

dim = size(mesh.elem, 2) - 1;
NT = size(mesh.elem, 1);

if ~isempty(codim);
	if (codim == 0)
		diam = simplex_diameter(mesh);
		return;
	end
	
	if isnumeric(codim)
		dim2 = dim - codim;
		toflip = false;
  elseif strcmpi(codim,'opposite')
		dim2 = dim - 1;
		toflip = true;
  else
    error('Unknown codimension');
  end
	
	% ASSEMBLE VECTOR MATRIX
	s = get_subsimplices(1:dim+1, dim2);
	ns = size(s,1);
	if dim2 == 0
		diam = ones(NT,ns);
		return;
	end
	edgelen = subsimplex_volume(mesh, dim - 1);
  
  
	edge = get_subsimplices(1:dim+1, 1);
	nedge = size(edge, 1);
	
	diam = zeros(NT, ns);
	for i=1:ns
    for j=1:nedge
      if any(s(i,:)==edge(j,1),2) && any(s(i,:)==edge(j,2),2)
        diam(:,i) = max(diam(:,i), edgelen(:,j));
      end
    end
	end
	
	if toflip, 
		diam = flip(diam,2); 
	end
elseif ~isempty(sub) 
	NS = size(sub,1);
	dim2 = size(sub,2)-1;
	if dim2 == 0
		diam = ones(NS,1);
		return;
  end
	
  diam = zeros(NS, 1);
  for i=1:dim2
    for j=i+1:dim2+1
      edgelen = sum(get_vec(mesh, sub(:,[i,j])).^2,2);
      diam = max(diam,edgelen);
    end
  end
  diam = sqrt(diam);
end

end