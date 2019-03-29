function mesh = micMeshG_X(x)
%MICMESHG_X Summary of this function goes here
%   Detailed explanation goes here

%% BDGMSH LOAD
[~, bdgmsh] = micBdmesh_X(x);

%% GMSH SETTINGS
options = struct;
options.lc = 0.15;
options.periodic = true;
options.optimize = true;

%% MESH GENERATION
mq = 0;
while mq < 0.3
	mesh = gmsh(bdgmsh, options);
	mq = min(simplex_quality(mesh));
	options.lc = options.lc / 1.2;
end
end