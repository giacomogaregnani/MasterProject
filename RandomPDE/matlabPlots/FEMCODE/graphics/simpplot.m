function simpplot(mesh, options)

%% DEPERIODIZE MESH
% does nothing if not periodic
figure('renderer','zbuffer')
[mesh, updates] = deperiodize(mesh);
mesh.bary = get_rc(mesh);

d=size(mesh.node,2);
dim = size(mesh.elem,2) -1;
NT = size(mesh.elem, 1);
N = size(mesh.node, 1);

if nargin<2, options = struct; end
if ~isfield(options,'nodedot'), options.nodedot = false; end
if ~isfield(options,'nodenum'), options.nodenum = false; end
if ~isfield(options,'elemnum'), options.elemnum = false; end
if ~isfield(options,'wire'), options.wire = false; end
if ~isfield(options,'facealpha'), options.facealpha = 1; end
if ~isfield(options,'edgealpha'), options.edgealpha = 1; end
if ~isfield(options,'fontsize'), options.fontsize = 10; end
if ~isfield(options,'view'), options.view = max(d,2); end
if ~isfield(options,'linewidth'), options.linewidth = 0.5; end
if ~isfield(options,'edgecolor'), options.edgecolor = 'black'; end
if ~isfield(options,'facecolor'), options.facecolor = 'white'; end
if ~isfield(options,'elevation'), options.elevation = 0; end

parsnodenum = {'fontname', 'times', ...
	'fontsize', options.fontsize, ...
	'color', 'black', ...
	'BackgroundColor',[.7 .9 .7], ...
	'FontWeight','bold', ...
	'Margin', 3, ...
	'VerticalAlignment', 'middle', ...
	'HorizontalAlignment','center', ...
	'EdgeColor','black', ...
	'LineWidth',1};
parsnodedot = {'linest','none', ...
	'marker', '.', ...
	'color', 'blue', ...
	'markersize', 24};
parselemnum={'fontname','times', ...
	'fontsize', options.fontsize, ...
	'FontWeight','bold', ...
	'VerticalAlignment', 'middle', ...
	'HorizontalAlignment','center'};
pars3dinface = {'edgecolor',options.edgecolor, ...
	'FaceColor','interp', ...
    'facealpha',options.facealpha};
pars3doutface= {'edgecolor',options.edgecolor, ...
	'FaceColor',options.facecolor, ...
    'facealpha',options.facealpha};
pars3dwire = {'edgecolor',options.edgecolor, ...
	'FaceColor','blue', ...
	'LineWidth', 2, ...
	'facealpha',options.facealpha, ...
	'edgealpha',options.edgealpha};

hold on;

if d == 1
	for it = 1:NT
		line(mesh.node(mesh.elem(it,:)), zeros(2,1));
	end
	if options.nodedot
		line(mesh.node, zeros(N,1), parsnodedot{:});
	end
	if options.nodenum
		for ip=1:N
			text(mesh.node(ip),0,num2str(updates(ip)),parsnodenum{:});
		end
	end
	if options.elemnum
		for it=1:NT
			text(mesh.bary(it), 0, num2str(it), parselemnum{:});
		end
	end
	axis off;
elseif d == 2
	trimesh(mesh.elem, mesh.node(:,1), mesh.node(:,2), ...
		repmat(options.elevation,size(mesh.node(:,1))), ...
		'facecolor', options.facecolor, ...
		'edgecolor',options.edgecolor, ...
		'linewidth', options.linewidth);
	if options.nodedot
		line(mesh.node(:,1), mesh.node(:,2), parsnodedot{:});
	end
	if options.nodenum
		for ip=1:N
			text(mesh.node(ip,1), mesh.node(ip,2), ...
				num2str(updates(ip)),parsnodenum{:});
		end
	end
	if options.elemnum
		for it=1:NT
			text(mesh.bary(it,1), mesh.bary(it,2), ...
			num2str(it), parselemnum{:});
		end
	end
	axis equal; axis off;
elseif (d == 3) && options.wire && (dim == 2)
	trisurf(mesh.elem, mesh.node(:,1),mesh.node(:,2),mesh.node(:,3), ...
		pars3dwire{:});
	axis equal; axis off;
elseif (d == 3) && options.wire && (dim == 3)
	face = auxstructure(mesh,'face');
	trisurf(face, mesh.node(:,1),mesh.node(:,2),mesh.node(:,3), ...
		pars3dwire{:});
	if options.nodenum
		for ip=1:N
			text(mesh.node(ip,1), mesh.node(ip,2), mesh.node(ip,3), ...
				num2str(updates(ip)),parsnodenum{:});
		end
	end
	if options.elemnum
		for it=1:NT
			text(mesh.bary(it,1), mesh.bary(it,2), mesh.bary(it,3), ...
			num2str(it), parselemnum{:});
		end
	end
	axis equal; axis off;
elseif (d == 3) && ~options.wire && (dim == 3)
	bdflag = auxstructure(mesh,'bdflag');
	if isfield(options,'expr')
		mesh.bary = get_rc(mesh);
		incl=eval(options.expr); % expr contains mesh.bary
		mesh.elem=mesh.elem(incl,:); % expr contains mesh.bary
		bdflag = bdflag(incl,:);
		if size(mesh.elem,1) > 0
			bdflag2 = auxstructure(mesh,'bdflag');
			bdflag2 = (bdflag2-bdflag)>0;
		else
			bdflag2 = false(0,d+1);
		end			
	else
		bdflag2 = false(NT, d+1);
	end
			
	leave = sum(bdflag,2) + sum(bdflag2,2) > 0;
	mesh.elem = mesh.elem(leave,:);
	bdflag  = bdflag (leave, :);
	bdflag2 = bdflag2(leave, :);
	[ind,~,mesh.elem(:)] = unique(mesh.elem(:));
	mesh.node = mesh.node(ind,:);
	
	if numel(mesh.elem)>0
		for i=1:4
			trisurf(mesh.elem(bdflag(:,i), setdiff(1:4,i)), ...
				mesh.node(:,1),mesh.node(:,2),mesh.node(:,3), ...
				pars3doutface{:});
			trisurf(mesh.elem(bdflag2(:,i), setdiff(1:4,i)), ...
				mesh.node(:,1),mesh.node(:,2),mesh.node(:,3), ...
				pars3dinface{:});
		end
	end
	axis equal; 
elseif (d == 3) && ~options.wire && (dim == 2)
	trisurf(mesh.elem, ...
		mesh.node(:,1), mesh.node(:,2), mesh.node(:,3), ...
		pars3dwire{:});
	trisurf(mesh.elem, ...
		mesh.node(:,1), mesh.node(:,2), mesh.node(:,3), ...
		pars3dwire{:});
	axis equal; axis off;
else
	error('Unimplemented dimension.');
end

view(options.view);
end
