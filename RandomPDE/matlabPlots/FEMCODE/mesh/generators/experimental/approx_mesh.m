function mesh = approx_mesh( bdmesh, options )

dim = size(bdmesh.node,2);

if nargin<2, options = struct; end
if ~isfield(options,'startN'), options.startN = 4; end
if ~isfield(options,'numiter'), options.numiter = 5; end
if ~isfield(options,'intersects'), options.intersects = struct; end
if ~isfield(options,'periodic'), options.periodic=false; end

%% FIRST MESH
mesh = structured_mesh(bdmesh.box, options.startN, ...
	struct('periodic', options.periodic));

%% ITERATIONS OF REFININT NEAR THE BOUNDARY
onbd = intersects(mesh, 'all', bdmesh, options.intersects);
for i=1:options.numiter
	[mesh,father] = bisect(mesh,onbd);
	onbd = onbd(father);
	onbd(onbd) = intersects(mesh, onbd, bdmesh, options.intersects);
end

%% CLEAR OUTSIDE ELEMENTS
bary = get_rc(mesh);
options.inside = 1;
[ ~, ~, inside] = dpoly(bary, bdmesh.node, bdmesh.elem, options);
mesh.elem = mesh.elem(inside,:);

%% CLEAR ELEMENTS WITH MORE THAN ONE BOUNDARY SIDES
mesh = cleanBdElem(mesh);

mesh = renumber(mesh); % CLEAR NODES THAT ARE NOT USED;
end

