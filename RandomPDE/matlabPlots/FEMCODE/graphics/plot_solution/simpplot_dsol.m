function simpplot_dsol(mesh,sol,options)
figure('renderer','zbuffer')
%   Copyright (C) 2004-2012 Per-Olof Persson. See COPYRIGHT.TXT for details.
[mesh, ~] = deperiodize(mesh);

if nargin<3
    options = struct;
end
if ~isfield(options,'edgecolor')
    options.edgecolor = 'none';
end
params = {'edgecolor',options.edgecolor};

d=size(mesh.node,2);
N = size(mesh.node,1);
NT = size(mesh.elem,1);
switch d
	case 1
		error('not yet implemented');
	case 2
		colormap(jet(256));
		trisurf(mesh.elem, mesh.node(:,1), mesh.node(:,2), ...
			sol(1:N),'FaceColor','flat','FaceVertexCData',sol(1:NT), ...
			params{:});
		view(2); colorbar; axis equal; axis off;
	case 3
		error('not yet implemented');
	otherwise
		error('Unimplemented dimension.');
end
