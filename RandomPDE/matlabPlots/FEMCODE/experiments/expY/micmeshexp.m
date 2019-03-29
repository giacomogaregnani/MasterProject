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
  mesh{i} = gmsh(bdgmsh{i},options);
end
mesh{6} = mesh{4};
mesh{6}.node = [mesh{6}.node; 1-mesh{6}.node];
N4 = size(mesh{4}.node,1);
mesh{6}.elem = [mesh{6}.elem; N4+mesh{6}.elem];
mesh{6} = fixorder(mesh{6});
for i=1:6
  mesh{i}.node = bsxfun(@plus,-[1/2,1/2],mesh{i}.node);
end


map = [2,5,5; 3,3,5; 4,3,2];
r   = [0,1,0; 0,2,1; 2,0,3];

fmesh = struct;
fmesh.node = zeros(0,2);
fmesh.elem = zeros(0,3);
N=0;
for i=1:3
  for j=1:3
    N = size(fmesh.node,1);
    fmesh.node = [fmesh.node; bsxfun(@plus, [j,4-i]-1/2, mesh{map(i,j)}.node*RC^(r(i,j)))];
    fmesh.elem = [fmesh.elem; N+mesh{map(i,j)}.elem];
  end
end


fmesh = remove_duplicate_nodes(fmesh);
fmesh.node = bsxfun(@plus, -[0.5, 0.5], 1/3*fmesh.node);
fmesh = periodize(fmesh,[-0.5,0.5,-0.5,0.5]);