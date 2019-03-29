function [varargout] = auxstructure(mesh, varargin)
% AUXSTRUCTURE provides auxiliary structures from mesh topology
%
% [varargout] = auxstructure(mesh, varargin) returns arrays describing the
% topological structure of the mesh requested by varargin.
%
% TOPOLOGICAL STRUCTURES
% 
% subsim = a lexicographically ordered list of all the subsimplices of
% codimension 1 (e.g. for 2D mesh it is a list of edges and for a 3D mesh
% it is a list of faces). The dimension of subsim is NE x (dim+1). This
% array defines a global index on subsimplices.
%
% elem2subsim = an array of dimension NT x (dim + 1). elem2subsim(i,j) is
% the index of the subsimplex that belongs to the element of global index i
% and is located opposite to the node with local index j.
%
% subsim2elem = an array of dimension NE x 4. Values in subsim2elem(i,1)
% and subsim2elem(i,2) are the indices of the two elements that share the
% sbusimplex with index i. If it is a boundary subsimplex, these values are
% the same. For j = 1,2, the subsimplex with index i is opposite the node
% with local number subsim2elem(i,2+j) in the element subsim2elem(i,j).
%
% neighbor = an arrat of dimension NT x (dim + 1). The value neighbor(i.j)
% is the index of the element that shares the subsimplex opposite the node
% with the local index j in the element number i.  If that subsimplex is a
% boundary subsimplex, neighbor(i,j) is set to i.
%
% bdsubsim = a (sub)list of boundary subsimplices
%
% isbdsubsim = a logical array of size NE x 1 that is true if the
% corresponding subsimplex is a boundary subsimplex.
%
% bdnode = a vector of nodes that are at the boundary
%
% is bdnode = a logical vector which is true for those nodes that are at
% the boundary
%
% bdflag = a logical array of size NT x (dim+1), where bdflag(i,j) is true
% if the subsimplex opposite to the node with local index j in the element
% number i is a boundary subsimplex.
%
% Further, for k = 1, ..., dim-1 one can ask the following outputs:
%
% simk (aliases: sim1 = edge, sim2 = face) = an array of dimension NE x 
% (k+1). For k=1 it is equivalent to 'subsim'.
% For a general k, the output is the list of all subsimplices of dimension
% k (therefore codimension dim-k) in lexicographical order. This array 
% defines a global index on subsimplices.
%
% elem2simk (aliases: elem2sim1 = elem2edge, elem2sim2 = elem2face) = an
% array of size NT x NS, where NS = binomial(dim+1,k+1) is the number of
% subsimplices of dimension k in a simplicial element of dimension dim. The
% row elem2simk(i,:) contains indices of all subsimplices (dimension k) of
% the element i. Their order is lexicographical in terms of local
% coordinates (Attention: For 2D, there is a difference between elem2subsim
% and elem2sim1, and that difference is in the order).
%
% bdsimk (aliases: bdsim1 = bdedge, bdsim2 = bdface) = only those rows
% from simk that represent boundary subsimplices.
%
% isbdsimk (aliases: isbdsim1 = isbdedge, isbdsim2 = isbdface) = a logical
% vector of size NE indicating boundary subsimplices.

%% INITIALIZATION
dim = size(mesh.elem,2)-1;
NT = size(mesh.elem,1);
N = size(mesh.node, 1);

varargout=cell(nargin-1,1);

maxo = 20+10*(dim-1) - 1;
wh   = zeros(maxo,1);
comp = false(maxo,1);

%% ALL INDICES
% 1  subsim       20 edge / sim1            40 sim3
% 2  elem2subsim  21 elem2edge / elem2sim1  41 elem2sim3
% 3  subsim2elem  22 bdedge / bdsim1        42 bdsim3
% 4  neighbor     23 isbdedge / isbdsim1    43 isbdsim3
% 5  bdsubsim     24 - 29 -                     ...
% 6  isbdsubsim   30 face / sim2               
% 7  bdnode       31 elem2face / elem2sim2
% 8  isbdnode     32 bdface / isbdface2
% 9  bdflag       33 isbdface / isbdsim2
% 10 - 19 -       34 - 39 -


%% READ INPUT
for i=1:nargin-1
	len = numel(varargin{i});
	if     strcmp(varargin{i},'subsim'),          c = 1; % OPPOSITE ORDER: 
	elseif strcmp(varargin{i},'elem2subsim'),     c = 2; % 1st SUBSIM IS 
	elseif strcmp(varargin{i},'subsim2elem'),     c = 3; % OPP. 1st NODE
	elseif strcmp(varargin{i},'neighbor'),        c = 4;
	elseif strcmp(varargin{i},'bdsubsim'),        c = 5;
	elseif strcmp(varargin{i},'isbdsubsim'),      c = 6;
	elseif strcmp(varargin{i},'bdnode'),          c = 7;
	elseif strcmp(varargin{i},'isbdnode'),        c = 8;
	elseif strcmp(varargin{i},'bdflag'),          c = 9;	
	elseif strcmp(varargin{i}(1:min(3,len)),'sim'),       % INCREASING
		c = 10*(str2double(varargin{i}(4:end))+1);        % ORDERING
	elseif strcmp(varargin{i}(1:min(8,len)),'elem2sim'),   
		c = 10*(str2double(varargin{i}(9:end))+1) + 1;
	elseif strcmp(varargin{i}(1:min(5,len)),'bdsim'),      
		c = 10*(str2double(varargin{i}(6:end))+1) + 2;
	elseif strcmp(varargin{i}(1:min(7,len)),'isbdsim'),    
		c = 10*(str2double(varargin{i}(8:end))+1) + 3;
	elseif (dim >= 2)	
		if strcmp(varargin{i},'edge'),            c = 20; 
		elseif strcmp(varargin{i},'elem2edge'),   c = 21; 
		elseif strcmp(varargin{i},'bdedge'),      c = 22;
		elseif strcmp(varargin{i},'isbdedge'),    c = 23;
		elseif (dim >= 3)
			if strcmp(varargin{i},'face'),            c = 30;
			elseif strcmp(varargin{i},'elem2face'),   c = 31;
			elseif strcmp(varargin{i},'bdface'),      c = 32;
			elseif strcmp(varargin{i},'isbdface'),    c = 33;
			else
				error('Unknown 3D or higher input');
			end
		else
			error('Unknown 2D input');
		end
	else
		error('unknown 1D input');
	end
	if (c <= maxo)
		comp(c) = true;
		wh(c)   = i;
		if (c >= 20) && ((mod(c,10) == 2) || (mod(c,10) == 3))
			comp(5) = true; % DEPENDENCIES
		end
	else
		error('Incompatible input');
	end
end

if sum(comp) == 0, return; end

if sum(comp(1:9)) > 0
	subsim = uint32(sort(get_subsimplices(mesh.elem,'opposite'),2));
	[subsim, i2, j] = unique(subsim,'rows','legacy');
end
if comp(1),
	varargout{wh(1)} = subsim; 
end
if sum(comp(2:9)) > 0
	i1(j((dim+1)*NT:-1:1)) = (dim+1)*NT:-1:1; 
	i1 = i1';
	if comp(2), 
		varargout{wh(2)} = uint32(reshape(j,NT,dim+1)); 
	end
	if sum(comp([3,4,9])) > 0
		k1 = ceil(i1/NT); t1 = i1 - NT*(k1-1);
		k2 = ceil(i2/NT); t2 = i2 - NT*(k2-1);
		if sum(comp([3,9])) > 0
			subsim2elem = uint32([t1,t2,k1,k2]);
			if comp(3),
				varargout{wh(3)} = subsim2elem;
			end
		end
		if sum(comp([4,9])) > 0
			idx = (i1 ~= i2);
			if comp(4)
				varargout{wh(4)} = ...
					uint32(accumarray([[t1(idx),k1(idx)]; ...
					[t2,k2]],[t2(idx);t1],[NT, dim+1]));
			end
			if comp(9)
				ind = find(i1 == i2);
				bdflag = accumarray([subsim2elem(ind,1), ...
					subsim2elem(ind,3)],true);
				if (size(bdflag,1)<NT) || (size(bdflag,2) < dim+1), 
					bdflag(NT,dim+1) = false; 
				end
				bdflag = bdflag>0;
				varargout{wh(9)} = bdflag;
			end
		end
	end
	if sum(comp(5:8)) > 0
		isbdsubsim = (i1 == i2);
		if comp(6), 
			varargout{wh(6)} = isbdsubsim; 
		end
		if sum(comp([5,7,8])) > 0
			bdsubsim = subsim(isbdsubsim,:);
			if comp(5) && (wh(5)>0), 
				varargout{wh(5)} = bdsubsim; 
			end
			if sum(comp([7,8])) > 0
				isbdnode = false(N,1);
				isbdnode(bdsubsim) = true;
				if comp(8), 
					varargout{wh(8)} = isbdnode; 
				end
				if comp(7) && (wh(7) > 0),
					varargout{wh(7)} = find(isbdnode); 
				end	
			end
		end
	end
end

for i=1:dim-1
	ind = (1+i)*10 : (1+i)*10 + 3;
	if sum(comp(ind(1:4))) > 0
		totalsimk = uint32(sort(get_subsimplices(mesh.elem,i),2));
		[simk, i2, j] = unique(totalsimk,'rows','legacy');
		if comp(ind(1)),
			varargout{wh(ind(1))} = simk;
		end
		if comp(ind(2))
			varargout{wh(ind(2))} = ...
				uint32(reshape(j,NT, nchoosek(dim+1,i+1)));
		end
		if sum(comp(ind(3:4))) > 0
			if i < dim-1
				bdsimk = ...
					unique(sort(get_subsimplices(bdsubsim,i),2),'rows');
				if comp(ind(3))
					varargout{wh(ind(3))} = bdsimk;
				end
				if comp(ind(4))
					[~, idx] = intersect(simk, bdsimk,'rows');
					isbdsimk = false(size(simk,1),1);
					isbdsimk(idx) = true;
					varargout{wh(ind(4))} = isbdsimk;
				end
			else
				clear i1;
				i1(j((dim+1)*NT:-1:1)) = (dim+1)*NT:-1:1;
				i1 = i1';
				isbdsimk = (i1 == i2);
				if comp(ind(3)) % bdsimk
					varargout{wh(ind(3))} = simk(isbdsimk,:);
				end
				if comp(ind(4)) % isbdsimk
					varargout{wh(ind(4))} = isbdsimk;
				end
			end	
		end
	end
end
end
