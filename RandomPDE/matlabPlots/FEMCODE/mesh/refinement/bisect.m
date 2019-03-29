function [mesh, father] = bisect(mesh, markedElem, father)
%% BISECT bisect an n-D triangulation.
% 
% [mesh] = bisect(mesh, markedElem) refines the current
% triangulation by bisecting marked elements and minimal neighboring
% elements to get a conforming and shape regular triangulation. Newest
% vertex bisection is implemented in 2D and longest edge bisection in 3D. 
% markedElem is a vector containing the indices of elements to be bisected. 
% It could be a logical vector of length size(elem,1) and even strin 'all'. 
% 
% - father(:, 1) is the reference of the element in which the refined
% elements lied before refining.
%
% TODO:
% 1. write the code for 1D

%% INPUT CONSTANTS
dim = size(mesh.elem,2) -1;
N = size(mesh.node,1); 
NT = size(mesh.elem,1);

%% INPUT READING
if (nargin == 1), markedElem = (1:NT)'; end
if isempty(markedElem), return; end
if strcmp(markedElem,'all'), markedElem = (1:NT)'; end
if islogical(markedElem), markedElem = find(markedElem); end
if nargin<3, father = (1:NT)'; end
	
if dim == 1
	error('NOT YET IMPLEMENTED - shame!');
elseif dim==2
	%% Construct auxiliary data structure
	[edge, elem2edge, neighbor] = ...
		auxstructure(mesh,'subsim','elem2subsim','neighbor');
	NE = size(edge,1);
	father(4*NT) = 0;
	
	%% Add new nodes
	isCutEdge = false(NE,1);
	while numel(markedElem)>0
		isCutEdge(elem2edge(markedElem,3)) = true;
		refineNeighbor = neighbor(markedElem,3);
		markedElem = ...
			refineNeighbor(~isCutEdge(elem2edge(refineNeighbor,3)));
	end
	edge2newNode = zeros(NE,1,'uint32');
	NN = sum(isCutEdge); % number of New Nodes
	edge2newNode(isCutEdge) = N+1:N+NN;
	mesh.node(N+1:N+NN,:) = get_rc(mesh,edge(isCutEdge,[1 2]));
	
	%% Refine marked elements
	Nb = 0; brother = zeros(3*NT,2,'uint32');
	for k = 1:2
		t = find(edge2newNode(elem2edge(:,3))>0);
		newNT = length(t);
		if (newNT == 0), break; end
		L = t; R = (NT+1:NT+newNT)';
		p1 = mesh.elem(t,1); p2 = mesh.elem(t,2); p3 = mesh.elem(t,3);
		p4 = edge2newNode(elem2edge(t,3));
		mesh.elem(L,:) = [p3, p1, p4];
		mesh.elem(R,:) = [p2, p3, p4];
		father(R) = father(t);
		brother(Nb+1:Nb+newNT,1) = L;
		brother(Nb+1:Nb+newNT,2) = R;
		elem2edge(L,3) = elem2edge(t,2);
		elem2edge(R,3) = elem2edge(t,1);
		
		% UPDATE BOUNDARY FLAG
		if isfield(mesh,'bdflag')
			subL = [0;3;1]; subR = [3;0;2];
			for j=1:numel(mesh.bdflag)
				[refind, refpos] = ismember(mesh.bdflag{j}(:,1),L);
				mesh.bdflag{j} = [mesh.bdflag{j}; 
					R(refpos(refind)),  subR(mesh.bdflag{j}(refind,2))];
				mesh.bdflag{j}(refind,2) = subL(mesh.bdflag{j}(refind,2));
				mesh.bdflag{j}(mesh.bdflag{j}(:,2) == 0,:) = [];
			end
		end
		
		NT = NT + newNT; Nb = Nb + newNT;
	end
	father = father(1:NT);
elseif (dim == 3)
	%% Pre-allocation
	mesh.node((2^dim+1)*N,1) = 0;
	mesh.elem(2^(dim-1)*NT,1) = 0;

	generation = zeros(N+6*NT,1,'uint8');

	father(2^(dim-1)*NT) = 0;
	
	%% Local Refinement
	%* Find new cutedges and new nodes
	%* Bisect all marked elements
	%* Find non-conforming elements
	%* Update generation of nodes
	
	cutEdge = zeros((2^dim)*N,3);      % cut edges
	nCut = 0;                          % number of cut edges
	nonConforming = true((2^dim)*N,1); % flag of non-conformity of edges
    
    mesh = label3(mesh);
	while ~isempty(markedElem)
		% switch element nodes such that 
		% elem(t,1:2) is the longest edge of t
		% OLD: [mesh,bdFlag] = label(mesh,markedElem,bdFlag);
		p = mesh.elem(markedElem,1:4);
		
		%% Find new cut edges and new nodes
		nMarked = length(markedElem); % number of marked elements
		p(nMarked,5) = 0;        % initialize new nodes
		if nCut == 0                  % In the first round, all marked
			idx = (1:nMarked)';       % elements introduce new cut edges;
		else                          % otherwise find existing cut edges
			ncEdge = find(nonConforming(1:nCut)); % all non-conform. edges
			nv2v = sparse([cutEdge(ncEdge,3);cutEdge(ncEdge,3)],...
				[cutEdge(ncEdge,1);cutEdge(ncEdge,2)],1,N,N);
			[i,j] = find(nv2v(:,p(:,1)).*nv2v(:,p(:,2)));
			p(j,5) = i;
			idx = find(p(:,5)==0);
		end
		if ~isempty(idx)                   % add new cut edges
			elemCutEdge = sort([p(idx,1) p(idx,2)],2); % all new cut edges
			% find(sparse) eliminates possible duplications in elemCutEdge
			[i,j] = find(sparse(elemCutEdge(:,1),elemCutEdge(:,2),1,N,N));
			nNew = length(i);              % number of new cut edges
			newCutEdge = nCut+1:nCut+nNew; % indices of new cut edges
			cutEdge(newCutEdge,1) = i;     % add cut edges
			cutEdge(newCutEdge,2) = j;     % cutEdge(:,1:2) two end nodes
			cutEdge(newCutEdge,3) = N+1:N+nNew; % cutEdge(:,3) middle point
			mesh.node(N+1:N+nNew,:) = get_rc(mesh,[i,j]); % add new nodes
			nCut = nCut + nNew;            % update number of cut edges
			N = N + nNew;                  % update number of nodes
			nv2v = sparse([cutEdge(newCutEdge,3);cutEdge(newCutEdge,3)],...
				[cutEdge(newCutEdge,1);cutEdge(newCutEdge,2)],1,N,N);
			[i,j] = find(nv2v(:,p(:,1)).*nv2v(:,p(:,2)));
			p(j,5) = i;
		end
		clear nv2v elemCutEdge
		
		%% Bisect marked elements
		idx = (generation(p(:,5)) == 0);
		if sum(idx > 0) == 1
			elemGeneration = ...
				max(transpose(generation(mesh.elem(markedElem(idx),:))));
		else
			elemGeneration = ...
				max(generation(mesh.elem(markedElem(idx),1:4)),[],2);
		end
		generation(p(idx,5)) = elemGeneration + 1;
        
        %% DEFINE NEW ELEMENTS AND RESOLVE LABELING
        % OLD: mesh.elem(markedElem,1:4) = [p4 p1 p3 p5];
		% OLD: mesh.elem(NT+1:NT+nMarked,1:4) = [p3 p2 p4 p5];
        me = false(NT,1); 
        me(markedElem) = true;
 
		% VERTICES OF CHILDERN 2
		[~, i2, i3] = indices((mesh.flag <= 3) | (mesh.flag == 10));
		mesh.elem(i2,1:4) = p(i3, [3,2,4,5]);
		[~, i2, i3] = indices(((mesh.flag >= 4) & (mesh.flag <= 6)) ...
			| (mesh.flag == 11));
		mesh.elem(i2,1:4) = p(i3, [2,4,3,5]);
		[~, i2, i3] = indices((mesh.flag >= 7) & (mesh.flag <= 9));
        mesh.elem(i2,1:4) = p(i3, [4,3,2,5]);
		
		% VERTICES OF CHILDREN 1
		[i1, ~, i3] = indices((mesh.flag == 1) | (mesh.flag == 4) | ...
			(mesh.flag == 7) | (mesh.flag == 10));
        mesh.elem(i1,1:4) = p(i3, [1,3,4,5]);
		[i1, ~, i3] = indices((mesh.flag == 2) | (mesh.flag == 5) | ...
			(mesh.flag == 8) | (mesh.flag == 11));
        mesh.elem(i1,1:4) = p(i3, [4,1,3,5]);
		[i1, ~, i3] = indices((mesh.flag == 3) | (mesh.flag == 6) | ...
			(mesh.flag == 9));
        mesh.elem(i1,1:4) = p(i3, [3,4,1,5]);
		
		% UPDATE LABELS AND FLAGS
		[ia1, ia2] = indices((mesh.flag == 1) | (mesh.flag == 5));
		[ib1, ib2] = indices(mesh.flag == 10);
		[ic1, ic2] = indices(mesh.flag == 11);
		
		% DIFFERENT LABELLING FOR Pf
		mesh.flag(markedElem) = 1;
		mesh.flag(NT+1:NT+nMarked) = 1;
		mesh.flag(ia1) = 10;
		mesh.flag(ia2) = 10;
		mesh.flag(ib1) = 4;
		mesh.flag(ib2) = 2;
		mesh.flag(ic1) = 2;
		mesh.flag(ic2) = 4;
		
        father(NT+1:NT+nMarked) = father(markedElem);
		NT = NT + nMarked;
		clear elemGeneration p
		
		%% Find non-conforming elements
		checkEdge = find(nonConforming(1:nCut)); % check non-conforming
		isCheckNode = false(N,1);                % edges
		isCheckNode(cutEdge(checkEdge,1)) = true; % check two end nodes of
		isCheckNode(cutEdge(checkEdge,2)) = true; % non-conforming edges
		isCheckElem = isCheckNode(mesh.elem(1:NT,1)) | ...
			isCheckNode(mesh.elem(1:NT,2)) | ...
			isCheckNode(mesh.elem(1:NT,3)) | ...
			isCheckNode(mesh.elem(1:NT,4));
		checkElem = find(isCheckElem); % all elements with checking nodes
		t2v = sparse(repmat(checkElem,4,1), ...
			mesh.elem(checkElem,:), 1, NT, N);
		[i,j] = find(t2v(:,cutEdge(checkEdge,1)) ...
			.* t2v(:,cutEdge(checkEdge,2)));
		markedElem = unique(i);
		nonConforming(checkEdge) = false;
		nonConforming(checkEdge(j)) = true;		
	end
	mesh.node = mesh.node(1:N,:);
	mesh.elem = mesh.elem(1:NT,:);
	father = father(1:NT);
else
	error('NOT IMPLEMENTED FOR dim > 3');
end

	function [i1, i2, i3] = indices(i0) 
		i1 = find(me & i0);
        i3 = find(i0(markedElem)); 
		i2 = NT + i3;
	end

mesh = cleanfields(mesh);
end

