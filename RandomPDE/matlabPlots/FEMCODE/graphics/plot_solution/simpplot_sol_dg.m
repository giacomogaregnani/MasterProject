function simpplot_sol_dg(mesh,sol,femspace,options)

d=size(mesh.node,2);

if nargin<4, options = struct; end
if ~isfield(options,'view'), options.view = 3; end

switch d
	case 2
            
            figure;
            vec1 = [1:3:3*size(mesh.elem,1);...
                    2:3:3*size(mesh.elem,1);...
                    3:3:3*size(mesh.elem,1)].';
            
            vec2 = [1:femspace.ldof:size(sol,1);...
                    2:femspace.ldof:size(sol,1);...
                    3:femspace.ldof:size(sol,1)];
            
            colormap(jet(256));
            trisurf(vec1, mesh.node(mesh.elem(:,:).',1), mesh.node(mesh.elem(:,:).',2), ...
                sol(vec2),'facecolor','interp','edgecolor','black');
            colorbar; axis equal; axis off;
            
	otherwise
		error('Unimplemented dimension.');
end
view(options.view);
