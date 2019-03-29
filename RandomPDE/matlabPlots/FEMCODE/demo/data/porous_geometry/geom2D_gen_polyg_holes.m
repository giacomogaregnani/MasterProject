function [bdmesh, bdgmsh]= geom2D_gen_polyg_holes(holes)
%GEOM2D_GEN_QUAD_HOLE Summary of this function goes here
%   Detailed explanation goes here

NH=numel(holes);

bdmesh.node = [-1/2, -1/2; 1/2, -1/2; 1/2, 1/2; -1/2, 1/2];
pos = 4;
h = zeros(NH,2);
for i=1:NH
	bdmesh.node= [bdmesh.node; holes{i}];
	h(i,1) = pos + 1;
	h(i,2) = pos + size(holes{i},1);
	pos = h(i,2);
end
bdmesh.elem=[1, 2; 2, 3; 3, 4; 4, 1];
for i=1:NH
	bdmesh.elem = [bdmesh.elem; [(h(i,1):h(i,2))', ...
		[(h(i,1)+1:h(i,2))'; h(i,1)]]];
end
bdmesh.box=[-0.5,0.5,-0.5,0.5];
bdmesh.per{1} = [1,3];
bdmesh.per{2} = [4,2];


bdmesh.bdnode{1,1} = 1; % lower left
bdmesh.bdnode{1,2} = 4; % upper left
bdmesh.bdnode{2,1} = 2; % lower right
bdmesh.bdnode{2,2} = 3; % upper right

bdmesh.bdnode{3,1} = [1, 2]; % down
bdmesh.bdnode{3,2} = [5, 3]; % up
bdmesh.bdnode{1,3} = [1, 4]; % left
bdmesh.bdnode{2,3} = [2, 3]; % right

bdmesh.bdelem{3,1} = 1;% down
bdmesh.bdelem{3,2} = 3; % up
bdmesh.bdelem{1,3} = 4; % left
bdmesh.bdelem{2,3} = 2; % right

bdgmsh.node = bdmesh.node;
bdgmsh.line = bdmesh.elem;
bdgmsh.periodicline = [1, -3; 2, -4];
bdgmsh.lineloop = {[1,2,3,4]};
for i=1:NH
	bdgmsh.lineloop{i+1} = h(i,1):h(i,2);
end
bdgmsh.planesurface = {1:NH+1};
bdgmsh.box=bdmesh.box;


end

