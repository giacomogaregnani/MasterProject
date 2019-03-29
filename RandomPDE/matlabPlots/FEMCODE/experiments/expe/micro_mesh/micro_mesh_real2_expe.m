function mesh2 = micro_mesh_real2_expe( p, vp )
%MICRO_MESH_REAL_EXPE Summary of this function goes here
%   Detailed explanation goes here

eps = vp.epsilon;
del = vp.delta;
R = del/eps;

kmin = ceil(p/eps - R/2 - 1);
kmax = floor(p/eps + R/2);

bdgmsh = struct;
bdgmsh.box = [-0.5,0.5,-0.5,0.5] * R;
bdgmsh.node = [-0.5,-0.5; 0.5,-0.5; 0.5,0.5; -0.5,0.5] * R;
bdgmsh.line = [1,2; 2,3; 3,4; 4,1];
if strcmp(vp.coupling,'periodic')
	bdgmsh.periodicline = [1, -3; 2, -4];
end
bdgmsh.lineloop{1} = [1,2,3,4];

bdmesh = struct;
bdmesh.node = zeros(0,2);
bdmesh.elem = zeros(0,2);
bdmesh2 = struct;
bdmesh2.node = zeros(0,2);
bdmesh2.elem = zeros(0,2);

for i = kmin(1):kmax(1)
	for j = kmin(2):kmax(2)
		x = eps * [i+1/2,j+1/2];
		[~,bdgmshloc] = pormat_expe(x,0);
		bdgmshloc.node = mdimop(bdgmshloc.node,[i+0.5,j+0.5]-p/eps,'+');
		if all(all(bdgmshloc.node(5:8,:) > bdgmsh.box(1)+0.1)) && ...
				all(all(bdgmshloc.node(5:8,:) < bdgmsh.box(2)-0.1))
			bdgmsh.node = [bdgmsh.node; bdgmshloc.node(5:8,:)];
			bdgmsh.line = [bdgmsh.line; bdgmsh.line(end-3:end,:)+4];
			bdgmsh.lineloop{end+1} = bdgmsh.lineloop{end} + 4;
			bdmesh2.node = [bdmesh2.node; bdgmshloc.node(5:8,:)];
			bdmesh2.elem = [bdmesh2.elem; [1,2; 2,3; 3,4; 4,1] + size(bdmesh2.elem,1)];
		end
		bdmesh.node = [bdmesh.node; bdgmshloc.node(5:8,:)];
		bdmesh.elem = [bdmesh.elem; [1,2; 2,3; 3,4; 4,1] + size(bdmesh.elem,1)];
	end
end
bdgmsh.planesurface{1} = 1:numel(bdgmsh.lineloop);



%% MESH GENERATION
options.lc = 0.3;
options.data_loc = '~/repos/experiments/auxdata/';
options.gmsh = '~/gmsh/bin/./gmsh';
options.periodic = strcmp(vp.coupling,'periodic');

%% MESH GENERATION
mq = 0;
while mq < 0.25
	mesh = gmsh(bdgmsh, options);
	mq = min(simplex_quality(mesh));
	options.lc = options.lc / 1.2;
end	

%% MESH REFINEMENT CLOSE TO THE BOUNDARY
int = dpoly(get_rc(mesh), bdmesh2.node, bdmesh2.elem) > 0.1;
for i=1:6
	onbd  = intersects(mesh, int, bdmesh);
	melem = int; melem(melem) = onbd; int=melem;
	[mesh,~,~,father] = bisect(mesh, melem);
	int = int(father);	
end
mesh2 = get_cleared_mesh(mesh,bdmesh);	

if strcmp(vp.coupling,'periodic')
	edge = auxstructure(mesh2,'edge');
	NT2 = size(mesh2.node,1);
	G = sparse( double(edge(:,1)), double(edge(:,2)), 1, NT2, NT2);
	G = G + G.'; %' make graph undirected
	[S, C] = graphconncomp( G ); % find connected components
	if S>1
		C = C';
		m = 0; num = 0; 
		for i=1:S
			newnum = sum(C==i);
			if newnum > num
				m=i;
				num = newnum;
			end
		end
		delnode = C ~= m;
		delelem = sum(delnode(mesh2.elem), 2) > 0;
		mesh2.elem(delelem,:) = [];
		[indi, ~, mesh2.elem(:)] = unique(mesh2.elem(:),'legacy');
		mesh2.node = mesh2.node(indi,:);
	end
end

%% ADD BOUNDARY CONDITIONS
tol = 10^(-5);
if strcmp(vp.coupling,'periodic')
	neufun = @(x)(false(size(x,1),1));
	dirfun = @(x)(true(size(x,1),1));
else
	neufun = @(x)(any(x > R/2-tol | x < -R/2 + tol,2));
	dirfun = @(x)( ~neufun(x));
end

mesh2.bdflag{1,1} = get_bdflag(mesh2, dirfun);
mesh2.bdflag{1,2} = mesh2.bdflag{1,1};
mesh2.bdflag{2,1} = get_bdflag(mesh2, neufun);
mesh2.bdflag{2,2} = mesh2.bdflag{2,1};

	function [mesh2, fat] = get_cleared_mesh(mesh2,bdmesh)
	bary = get_rc(mesh2);
	opt.inside = 1;
	[ ~, ~, inside] = dpoly(bary, bdmesh.node, bdmesh.elem, opt);
	mesh2.elem = mesh2.elem(inside,:);
	fat = find(inside);
	
	%% CLEAR ELEMENTS WITH MORE THAN ONE BOUNDARY SIDES
	NT = size(mesh2.elem);
	NTold = size(mesh2.elem)+1;
	while sum(NT<NTold)
		[isbdsubsim,elem2subsim] = ...
			auxstructure(mesh2,'isbdsubsim','elem2subsim');
		idx = (sum(isbdsubsim(elem2subsim),2) <= 1);
		mesh2.elem = mesh2.elem(idx,:);
		fat = fat(idx);
		NTold = NT;
		NT = size(mesh2.elem);
	end
	
	%% CLEAR NODES THAT ARE NOT USED;
	[ind, ~, mesh2.elem(:)] = unique(mesh2.elem(:),'legacy');
	mesh2.node = mesh2.node(ind,:);
	end

end

