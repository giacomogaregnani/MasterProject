load('prod_mesh.mat');
mesh = struct('node',p,'elem',t);
bdmesh = struct('node', mesh.node, 'elem',  auxstructure(mesh,'bdface'));
bdmesh.elem(any(abs(get_rc(bdmesh))<0.001,2) | any(abs(get_rc(bdmesh)-1)<0.001,2),:) = []; 
bdmesh = renumber(bdmesh);
simpplot(bdmesh);

N = size(bdmesh.node,1);
s = N+1+div*(0:2);
z = s(3)+div+(0:3);

bdgmsh = struct;
bdgmsh.node = [bdmesh.node;
  bsxfun(@plus, bdmesh.node(1:div,:), [2,0,0]);
  bsxfun(@plus, bdmesh.node(div+1:2*div,:), [0,2,0]);
  bsxfun(@plus, bdmesh.node(2*div+1:3*div,:), [0,0,2]);
  3*eye(3)];

