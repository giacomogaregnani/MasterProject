function tricontour2(mesh, sol, M, lw, elev)

% Contouring for functions defined on triangular meshes
%
%   TRICONTOUR(p,t,F,N)
%
% Draws contours of the surface F, where F is defined on the triangulation
% [p,t]. These inputs define the xy co-ordinates of the points and their
% connectivity:
%
%   P = [x1,y1; x2,y2; etc],            - xy co-ordinates of nodes in the 
%                                         triangulation
%   T = [n11,n12,n13; n21,n23,n23; etc] - node numbers in each triangle
%
% The last input N defines the contouring levels. There are several
% options:
%
%   N scalar - N number of equally spaced contours will be drawn
%   N vector - Draws contours at the levels specified in N
%
% A special call with a two element N where both elements are equal draws a
% single contour at that level.
%
%   [C,H] = TRICONTOUR(...)
%
% This syntax can be used to pass the contour matrix C and the contour
% handels H to clabel by adding clabel(c,h) or clabel(c) after the call to
% TRICONTOUR.
%
% TRICONTOUR can also return 3D contours similar to CONTOUR3 by adding
% view(3) after the call to TRICONTOUR.
%
% Type "contourdemo" for some examples.
%
% See also, CONTOUR, CLABEL

% This function does NOT interpolate back onto a Cartesian grid, but
% instead uses the triangulation directly.
%
% If your going to use this inside a loop with the same [p,t] a good
% modification is to make the connectivity "mkcon" once outside the loop
% because usually about 50% of the time is spent in "mkcon".
%
% Darren Engwirda - 2005 (d_engwirda@hotmail.com)
% Updated 15/05/2006

[mesh, updates] = deperiodize(mesh);
sol = sol(updates);

d=size(mesh.node,2);
dim = size(mesh.elem,2) - 1;
%N = size(mesh.node,1);
%M = 20;

if nargin<4
	lw = 0.5;
end

if nargin<5
  elev = 2;
end

if numel(M) > 1, 
	v = M;
else
	v = linspace(min(sol),max(sol),M);
	v([1,M]) = [];
end

[subsim, elem2subsim] = ...
	auxstructure(mesh,'subsim', 'elem2subsim');
NE=size(subsim,1);
N=size(mesh.node,1);
for i=1:numel(v)
	val = v(i);
	cutedge = (prod((sol(subsim) - val),2) < 0);
	cutpoints = zeros(NE+N+N,d);
	cutpoints(NE+1:end,:) = [mesh.node; mesh.node];
	ve1 = sol(subsim(cutedge,1));
	ve2 = sol(subsim(cutedge,2));
	for k=1:d
		cutpoints(cutedge,k) = (...
			mesh.node(subsim(cutedge,1),k) .* (ve2 - val) + ...
			mesh.node(subsim(cutedge,2),k) .* (val - ve1)) ...
			./ (ve2-ve1);
	end

	A = sparse([],[],[],NE+N+N,NE+N+N);
	tr = [1,2,3; 1,3,2; 2,3,1];
	scutelem = sum(cutedge(elem2subsim),2) == 1;
	for j=1:3
		cutelem = cutedge(elem2subsim(:,tr(j,1))) & ...
			cutedge(elem2subsim(:,tr(j,2)));
		B = sparse(double(elem2subsim(cutelem,tr(j,1))), ...
			double(elem2subsim(cutelem,tr(j,2))),1,NE+2*N,NE+2*N);
		A = A + B + B';
		cutelem = scutelem & cutedge(elem2subsim(:,j));
		B = sparse(double(elem2subsim(cutelem,j)), ...
			double(NE + mesh.elem(cutelem,j)), ...
			1,NE+2*N,NE+2*N);	
		A = A + B + B';
	end
	fulle = sum(sol(subsim) == val,2) == 2;
	B = sparse(NE + N + double(subsim(fulle,1)), ...
		NE + N + double(subsim(fulle,2)), ...
		1, NE+2*N, NE+2*N);
	A = A + B + B';
	
	[cutedgelist, ~] = find(A);
	cutedgelist = unique(cutedgelist)';
	cutedge = false(NE+2*N,1);
	cutedge(cutedgelist) = true;
	for j=cutedgelist
		if cutedge(j)
			nbrs = find(A(:,j));
			edgelist = j;
			newedge = nbrs(1);
			prevedge = j;
			while newedge > 0 && newedge ~= j 
				edgelist = [edgelist; newedge];
				newnbrs = find(A(:,newedge));
				if numel(newnbrs) == 1
					newedge = 0;
				else
					backup = newedge;
					newedge = setdiff(newnbrs, prevedge);
					prevedge = backup;
				end
			end
			if newedge == 0 && numel(nbrs) > 1
				newedge = nbrs(2);
				prevedge = j;
				while newedge > 0 && newedge ~= j
					edgelist = [newedge; edgelist];
					newnbrs = find(A(:,newedge));
					if numel(newnbrs) == 1
						newedge = 0;
					else
						backup = newedge;
						newedge = setdiff(newnbrs, prevedge);
						prevedge = backup;
					end
				end
			end
			cutedge(edgelist) = false;
			if newedge == 0
				if d==2
					line(cutpoints(edgelist,1), ...
						cutpoints(edgelist,2), ...
						repmat(elev,size(edgelist)),...
						'color','black','linewidth',lw)
				elseif d==3
					line(cutpoints(edgelist,1), ...
						cutpoints(edgelist,2), ...
						cutpoints(edgelist,3), ...
						'color','black','linewidth',lw)
				end
			else
				if d==2
					patch(cutpoints(edgelist,1), ...
						cutpoints(edgelist,2),...
						repmat(elev,size(edgelist)), ...
						[1,1,1], ...
						'linewidth',lw,'facecolor','none')
				elseif d==3
					patch(cutpoints(edgelist,1), ...
						cutpoints(edgelist,2),...
						cutpoints(edgelist,3), ...
						[1,1,1], ...
						'linewidth',lw,'facecolor','none')
				end
			end
		end
	end
end

