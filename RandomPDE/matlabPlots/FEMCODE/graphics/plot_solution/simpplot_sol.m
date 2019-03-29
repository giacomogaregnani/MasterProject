function simpplot_sol(mesh,f,options)
figure('renderer','zbuffer')
[mesh, updates] = deperiodize(mesh);
if isnumeric(f)
	f = f(updates);
else
	f = f(mesh.node);
end

d=size(mesh.node,2);
dim=size(mesh.elem,2)-1;
N = size(mesh.node,1);

if nargin<3, options = struct; end
if ~isfield(options,'view'), options.view = max(d,2); end
if ~isfield(options,'colorbar'), options.colorbar = true; end
if ~isfield(options,'equal'), options.equal = true; end
if ~isfield(options,'off'), options.off = true; end
if ~isfield(options,'edgecolor'), options.edgecolor = 'interp'; end

hold on;

switch d
	case 1
		[sorted, a] = sort(mesh.node);
		plot(sorted, f(a));
	case 2
		colormap(jet(256));
		trisurf(mesh.elem, mesh.node(:,1), mesh.node(:,2), ...
			f(1:N),'facecolor','interp','edgecolor',options.edgecolor);
	case 3
		if nargin>2 && isfield(options,'expr')
			mesh.bary = get_rc(mesh);
			incl=eval(options.expr); % expr contains mesh.bary
			mesh.elem=mesh.elem(incl,:); % expr contains mesh.bary
		end
		if dim == 3
			bdflag = auxstructure(mesh,'bdflag');
			leave = sum(bdflag,2)>0;
			mesh.elem = mesh.elem(leave,:);
			[ind,~,mesh.elem(:)] = unique(mesh.elem(:));
			mesh.node = mesh.node(ind,:);
			f = f(ind);
			
			bdflag = auxstructure(mesh,'bdflag');
			colormap(jet(1024));
			for i=1:4
				if numel(mesh.elem)>0
					h=trisurf(mesh.elem(bdflag(:,i), setdiff(1:4,i)), ...
						mesh.node(:,1),mesh.node(:,2),mesh.node(:,3));
					set(h,'edgecolor',options.edgecolor,'FaceColor','interp', ...
						'FaceVertexCData',f);
				end
				hold on
			end
		elseif dim == 2
			[ind,~,mesh.elem(:)] = unique(mesh.elem(:));
			mesh.node = mesh.node(ind,:);
			f = f(ind);
			h=trisurf(mesh.elem, ...
				mesh.node(:,1),mesh.node(:,2),mesh.node(:,3));
			set(h,'edgecolor',options.edgecolor,'FaceColor','interp', ...
				'FaceVertexCData',f);
		end
	otherwise
		error('Unimplemented dimension.');
end
if options.colorbar, colorbar; end
if options.off, axis off; end
if options.equal, axis equal; end
view(options.view);
