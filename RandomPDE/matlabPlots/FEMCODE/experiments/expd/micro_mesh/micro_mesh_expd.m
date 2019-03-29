function mesh = micro_mesh_expd( p, varargin )
%MICRO_MESH_EXPA Summary of this function goes here
%   Detailed explanation goes here

%% BDGMSH LOAD
[~,bdgmsh] = pormat_expd(p,0);

%% GMSH SETTINGS
options.lc = 0.35;
options.data_loc = '~/repos/experiments/auxdata/';
options.gmsh = '~/gmsh/bin/./gmsh';

%% MESH GENERATION
mesh = gmsh(bdgmsh, options);

end

