function mesh = mesh_polyhedron(bdmesh, options, pfix, pstart )
%MESH_POLYHEDRON meshes a polyhedron using DISTMESH algorithm
%
% mesh = mesh(bdmesh) produces a mesh in the polygon described by the
% boundary bdmesh. 

if nargin < 2, options = struct; end
if nargin < 3, pfix = []; end
if nargin < 4, pstart = []; end
if ~isfield(options,'h'), options.h = 0.1; end
if ~isfield(options,'maxit'), options.maxit = 20; end
if ~isfield(options,'periodic'), options.periodic = false; end
if ~isfield(options,'fix2'), options.fix2 = false; end
if ~isfield(options,'fix1'), options.fix1 = false; end

d=size(bdmesh.node,2);
dim = size(bdmesh.elem,2)-1;

if (dim+1 ~= d)
	error('boundary does hot have the proper dimension');
end

h = options.h;
box=[min(bdmesh.node); max(bdmesh.node)];

%% FIX VERTICES (2D) or EDGES (3D)
if options.fix2
	if (d==2)
		pfix = bdmesh.node;
	elseif (d==3)
		pfix = bdmesh.node;
		edge = auxstructure(bdmesh,'sim1');
		for i=1:size(edge,1)
			% force nodes between edge(i,:) and edge(i+1,:)
			n1 = bdmesh.node(edge(i,1),:);
			n2 = bdmesh.node(edge(i,2),:);
			dist = sqrt(sum((n2-n1).^2));
			numstep = ceil(dist/h/1.6);
			step = (n2-n1)/numstep;
			pfix = [pfix; repmat(n1,[numstep-1,1]) + ...
				repmat((1:numstep-1)',[1,3]) .* ...
				repmat(step, [numstep-1,1])];
		end
	else
		error('FIX 2 not implemented for this dimension');
	end
end

%% FIX EDGES (2D) or FACES (3D)
if options.fix1
	if (d == 2)
		if ~options.periodic
			pfix = bdmesh.node;
			for i=1:size(bdmesh.elem,1)
				% force nodes between bdmesh.node(i,:) and bdmesh.node(i+1,:)
				n1 = bdmesh.node(bdmesh.elem(i,1),:);
				n2 = bdmesh.node(bdmesh.elem(i,2),:);
				dist = sqrt(sum((n2-n1).^2));
				numstep = ceil(dist/h);
				step = (n2-n1)/numstep;
				pfix = [pfix; repmat(n1,[numstep-1,1]) + ...
					repmat((1:numstep-1)',[1,2]) .* ...
					repmat(step, [numstep-1,1])];
			end
		end
	else
		error('FIX 2 not implemented for this dimension');
	end

end


%ptol=.001; deps=sqrt(eps)*h;
ttol=.1; L0mult=1+.4/2^(d-1); deltat=.1; geps=0.1*h; 
densityctrlfreq=10;

% 1. Create initial distribution in bounding box
if isempty(pstart)
	if d==1
		p=(box(1):h:box(2))';
	else
		cbox=cell(1,d);
		for ii=1:d
			cbox{ii}=box(1,ii):h:box(2,ii);
		end
		pp=cell(1,d);
		[pp{:}]=ndgrid(cbox{:});
		p=zeros(numel(pp{1}),d);
		for ii=1:d
			p(:,ii)=pp{ii}(:);
		end
	end
else
	p = pstart;
end

% 2. Remove points outside the region, apply the rejection method
options.inside = 1;
[distance,~,inside,~,~] = dpoly(p, bdmesh.node, bdmesh.elem, options);
p=p((inside .* (distance > (0.3*h)^2)) > 0, :);
if ~isempty(pfix)
	p=[pfix; p];
end
N=size(p,1);

count=0;
count2=1;
p0=inf;
keep = false;
while (count2<options.maxit) || keep % 3. Retriangulation by Delaunay
	keep = false;
	if max(sqrt(sum((p-p0).^2,2)))>ttol*h
		p0=p;
		t=delaunayn(p);
		pmid=zeros(size(t,1),d);
		for ii=1:d+1
			pmid=pmid+p(t(:,ii),:)/(d+1);
		end
		[distance,~,inside,~,~] = ...
			dpoly(pmid, bdmesh.node, bdmesh.elem, options);
		t=t((inside + (distance > geps^2)) > 1, :);
		% 4. Describe each edge by a unique pair of nodes
		pair=zeros(0,2);
		localpairs=nchoosek(1:d+1,2);
		for ii=1:size(localpairs,1)
			pair=[pair;t(:,localpairs(ii,:))];
		end
		pair=unique(sort(pair,2),'rows');
		count=count+1;
		keep = true;
	end
		% 6. Move mesh points based on edge lengths L and forces F
		bars=p(pair(:,1),:)-p(pair(:,2),:);
		L=sqrt(sum(bars.^2,2));
	
		L0=h*L0mult*(sum(L.^d)/(size(L,1)*h^d))^(1/d);
	
		F=max(L0-L,0);
		Fbar=[bars,-bars].*repmat(F./L,1,2*d);
		dp=full(sparse(pair(:,[ones(1,d),2*ones(1,d)]), ...
			ones(size(pair,1),1)*[1:d,1:d], ...
			Fbar,N,d));
		if ~isempty(pfix)
			dp(1:size(pfix,1),:)=0;
		end
		p=p+deltat*dp;

	% Density control - remove points that are too close
	if mod(count2,densityctrlfreq)==0 && any(L0>2*L)
		%p(reshape(pair(L0>2*L,1),[],1),:)=[]; before 
		pair = sort(pair,2);
		p(reshape(pair(L0>2*L,2),[],1),:)=[]; 
		N=size(p,1); p0=inf;
		keep = true;
	end
	
	% 7. Bring outside points back to the boundary
	[ ~, closest, inside, ~, ~ ] = ...
		dpoly( p, bdmesh.node, bdmesh.elem, options);
	inside(1:size(pfix,1)) = false;
	p(~inside,:) = closest(~inside,:);
	
	% 8. Termination criterion
	%maxdp=max(deltat*sqrt(sum(dp(d<-geps,:).^2,2)));
	%if maxdp<ptol*h, break; end
	count2= count2+1;
end

mesh.elem = t;
mesh.node = p;
mesh.periodic = false;
end
