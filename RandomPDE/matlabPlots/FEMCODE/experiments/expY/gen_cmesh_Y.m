function mesh = gen_cmesh_Y

options = struct;
options.lc = 0.1;
options.periodic = false;
options.optimize = 1;

ang = linspace(0,pi/2,9)';

clear bdgmsh mesh
bdgmsh{1}.node = [0,0; 1,0; 1,1; 0,1];
bdgmsh{2}.node = [0,0; 1,0; 1,0.5; 0,0.5];
bdgmsh{3}.node = [1/2*cos(ang), 1/2*sin(ang); 0,1; 1,1; 1,0];
bdgmsh{4}.node = [1/2*cos(ang), 1/2*sin(ang); 0,0];
bdgmsh{5}.node = [1/2*cos(ang), 1/2*sin(ang); 0,1; 1-1/2*cos(ang), 1-1/2*sin(ang);  1,0];

RC = [0,1; -1,0];

for i=1:5
  N = size(bdgmsh{i}.node,1);
  bdgmsh{i}.line = [(1:N)', [(2:N)'; 1]];
  bdgmsh{i}.lineloop = {(1:N)};
  bdgmsh{i}.planesurface = {1};
  meshes{i} = gmsh(bdgmsh{i},options);
end
meshes{6} = meshes{4};
meshes{6}.node = [meshes{6}.node; 1-meshes{6}.node];
N4 = size(meshes{4}.node,1);
meshes{6}.elem = [meshes{6}.elem; N4+meshes{6}.elem];
meshes{6} = fixorder(meshes{6});
for i=1:6
  meshes{i}.node = bsxfun(@plus,-[1/2,1/2],meshes{i}.node);
end


map = [2,5,5; 3,3,5; 4,3,2];
r   = [0,1,0; 0,2,1; 2,0,3];

mesh = struct;
mesh.node = zeros(0,2);
mesh.elem = zeros(0,3);
N=0;
for i=1:3
  for j=1:3
    N = size(mesh.node,1);
    mesh.node = [mesh.node; bsxfun(@plus, [j,4-i]-1/2, meshes{map(i,j)}.node*RC^(r(i,j)))];
    mesh.elem = [mesh.elem; N+meshes{map(i,j)}.elem];
  end
end


mesh = remove_duplicate_nodes(mesh);
mesh.node = bsxfun(@plus, -[0.5, 0.5], 1/3*mesh.node);
mesh = periodize(mesh,[-0.5,0.5,-0.5,0.5]);

end

