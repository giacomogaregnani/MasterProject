function mesh = micMesh_e( p, varargin )
%MICRO_MESH_EXPA Summary of this function goes here
%   Detailed explanation goes here

%% BDGMSH LOAD
[~,bdgmsh] = pormat_expe(p,0);

%% GMSH SETTINGS

options.lc = 0.3;
options.data_loc = '~/repos/experiments/auxdata/';
options.gmsh = '~/gmsh/bin/./gmsh';

%% MESH GENERATION
mq = 0;
while mq < 0.6
	mesh = gmsh(bdgmsh, options);
	mq = min(simplex_quality(mesh));
	options.lc = options.lc / 1.2;
end	
end

