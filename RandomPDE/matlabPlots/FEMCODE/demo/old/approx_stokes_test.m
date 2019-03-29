function approx_stokes_test
%% set porous geometry parameters
%clear fg sector options;
%fg.epsilon = 0.1;
%fg.handle = @cg2d_1;
%sector.centre = [0.47,0.22];
%sector.delta = fg.epsilon;

%% generate boundary
r = 0.1;
bdmesh = par2geom_expe(r);
load('ahom_expe.mat','cs');
itensor = [ppval(cs{1,1},r), ppval(cs{1,2},r); ppval(cs{2,1},r), ppval(cs{2,2},r)];

%% set approx geometry parameters
clear options;
options.startN = 8;
options.numiter = 0;
options.int_approx = false;
options.periodic = true;

vpm = struct;
vpm.f = @fstokesmicro1;
vpm.elemtype = 'p2';
vpm.pelemtype = 'p1';
vpm.a = 1;
vpm.bc = 'zero_dirichlet';

%% generate approx. geometry

mesh = structured_mesh(bdmesh.box, options.startN, ...
	struct('periodic', options.periodic));

k=1; dim = size(mesh.elem,2) - 1;
while true
	%% CLEAN GEOMETRY
	[mesh2, fat] = get_cleared_mesh(mesh);
    
	%% SOLVE
	[usol, ufemspace, psol, pfemspace, ...
        mesh2, vpm] = stokes(mesh2, vpm); % no options

    average(k,:) = get_integral(mesh2, ufemspace, usol);
	ndof(k) = dim*ufemspace.ndof + pfemspace.ndof;
	nelem(k) = size(mesh2.elem,1);
	nnode(k) = size(mesh2.node,1);
	
    %% ESTIMATE
	res = stokes_residual(mesh2, vpm, ufemspace, usol, pfemspace, psol);
	abserr(k) = sum(res, 1);	
	fprintf(1,'iter %d, error %f\n',k,abserr(k));
	
	if abserr(k) < 0.0001, % STOPPING CRITERIA
		break;
	end
	
	markedElements = markelem(res); % MARK
	markedElements = fat(markedElements);
	mesh = bisect(mesh, markedElements); % REFINE
	k=k+1;  % ADVANCE
end



%% plot
figure;
simpplot(mesh2);
simpplot(bdmesh);

	function [mesh2, fat] = get_cleared_mesh(mesh2)
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
	[i, ~, mesh2.elem(:)] = unique(mesh2.elem(:),'legacy');
	mesh2.node = mesh2.node(i,:);
	end
end