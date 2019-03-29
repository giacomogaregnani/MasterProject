function mesh = micro_mesh_expc( p, varargin )
%MICRO_MESH_EXPA Summary of this function goes here
%   Detailed explanation goes here

%% BDGMSH LOAD
[~,bdgmsh] = pormat_expc(p,0);

%% GMSH SETTINGS
options.lc = 0.2;
options.data_loc = '~/repos/experiments/auxdata/';
options.gmsh = '~/gmsh/bin/./gmsh';

%% MESH GENERATION
mq = 0;
while mq < 0.3
	mesh = gmsh(bdgmsh, options);
	mq = min(simplex_quality(mesh));
	options.lc = options.lc / 1.2;
end	
end

