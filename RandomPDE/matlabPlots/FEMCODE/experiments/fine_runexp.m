%% USER INPUT

% required - normally changed
expName  = 'X';
umicElem = 'p2';
pmicElem = 'p1';
M = 16; % refinement level
lc = 0.05;
saveDir = '~/repos/experiments/comp/';
map = [true(3,4); false(1,4)];


%% DO NOT TOUCH
mesh = composite_mesh(@micBdmesh_X,map,M,struct('periodic',true,'lc',lc));
bary = get_rc(mesh);

nodes = [0,0; 1,0; 2,0; 2,1; 3,2; 3,3; 2,3; 1,3; 2,4; 1,4; 0,4; 0,3; 1,2; 2,2; 1,1; 0,1];
edges = [(1:16)', [(2:16)'; 1]];

[~, ~, inside, ~, ~] = dpoly(bary, nodes, edges, struct('inside',true));


mesh.elem = mesh.elem(inside,:); % throw away outside elements
mesh = cleanBdElem(mesh); % REMOVE ELEMENTS WITH MORE THAN ONE BOUNDARY SIDES
mesh = renumber(mesh); % remove unused nodes


saveFile = [saveDir expName '/mesh' num2str(M) '.mat'];
save(saveFile,'mesh');