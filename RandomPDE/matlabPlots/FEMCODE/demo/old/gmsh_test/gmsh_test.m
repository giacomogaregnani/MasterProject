%options.data_loc = 'D:/repos/fehmm/experiments/auxdata/';
%options.gmsh = '"D:/Program Files/gmsh/gmsh.exe"';
options.lc = 0.2;
[bdmesh, bdgmsh] = geom2D_exp([0.05,0.4,0.05,0.4],1/8);
mesh = gmsh(bdgmsh, options);
mesh = label(mesh);
%simpplot(mesh, struct('facealpha',0.3));
simpplot(mesh);
