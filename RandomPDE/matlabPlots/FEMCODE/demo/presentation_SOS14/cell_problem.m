function cell_problem

%% Define coarse problem
vp = struct('elemtype', 'p2', 'pelemtype', 'p1', 'bc', [0,0], 'a', 1, 'f', @(x)(repmat([0,-1], size(x,1),1)));
point = [0,0]; bdmesh = pormat_expc(point,0); mesh = micro_mesh_expc(point); mesh.bdflag = 'dirichlet';

cmax = 0.7;  cmin = -0.6;
mag = '-m8'; aa = '-a2'; bgcolor = [0.1,0.3,1];
scale = 3;

N = [1];
for i=1:numel(N)
	if i<2, lastN=0; 
	else lastN=N(i-1);
	end
	for j=lastN+1:N(i)
		mesh = uniformrefine(mesh);
	end
		
	%% plot mesh
	if N(i) < 2
		init_figure;
		simpplot(mesh, struct('facecolor',bgcolor));
		export_fig(['mesh' num2str(N(i)) '.png'],'-transparent',mag,aa);
		close
	end
	
	%% solve coarse problem
	[usol, ~, psol, ~, mesh, vp] = stokes(mesh, vp);
	psol(psol > cmax) = cmax; psol(psol < cmin) = cmin;
	
	%% plot coarse pressure
	init_figure;
	simpplot_sol(mesh,psol); 
	if N(i)<2, hold on; simpplot(mesh, struct('elevation',1,'facecolor','none')); end
	h = colorbar; pos = get(h,'position'); pos(1) = 0.80; pos(2) = 0.2; pos(4) = 0.6; set(h,'position',pos); caxis([min(psol),max(psol)]); colormap(CMRmap(256));
	export_fig(['pressure' num2str(N(i)) '.png'],'-transparent',mag,aa);
	close
	
	%% plot coarse velocity
	if N(i)<2
		[mesh2,update] = deperiodize(mesh);
		usol2 = usol(update,:);
		
		init_figure;
		quiver(mesh2.node(:,1), mesh2.node(:,2), scale*usol2(:,1), scale*usol2(:,2), 'color','white','autoscale','off');		
		hold on; simpplot(mesh, struct('elevation',-1,'facecolor',bgcolor));		
		axisset;
		export_fig(['velocity' num2str(N(i)) '.png'],'-transparent',mag,aa);
		close
		
		% "fine velocity"
		if N(i) == 1
			init_figure;
			quiver(mesh2.node(:,1), mesh2.node(:,2), scale*usol2(:,1), scale*usol2(:,2), 'color','white','autoscale','off');
			hold on; 
			simpplot(mesh, struct('elevation',-2,'edgecolor',bgcolor ,'facecolor',bgcolor));
			simpplot(bdmesh, struct('elevation',-1,'facecolor','none'));			
			axisset;
			export_fig('velocity_fine.png','-transparent',mag,aa);
			close
			
			init_figure;
			bdnode = auxstructure(mesh2,'bdnode');
			quiver(mesh2.node(bdnode,1), mesh2.node(bdnode,2), scale*usol2(bdnode,1), scale*usol2(bdnode,2), 'color','white','autoscale','off');
			simpplot(mesh, struct('elevation',-2,'edgecolor',bgcolor ,'facecolor',bgcolor));
			hold on; simpplot(bdmesh, struct('elevation',-1,'facecolor','none'));	
			axisset;
			export_fig('problem_setting.png','-transparent',mag,aa);
			close
			
			init_figure;
			simpplot(mesh, struct('elevation',-2,'edgecolor',bgcolor ,'facecolor',bgcolor));
			hold on; simpplot(bdmesh, struct('elevation',-1,'facecolor','none'));	
			axisset;
			export_fig('bdmesh.png','-transparent',mag,aa);
			close
		end
	end
end
	function init_figure
	figure; set(gcf,'color','black','position',[100,100,270,180]);
	end

	function axisset
		xlim([-0.6,0.6]); ylim([-0.6,0.6]); axis equal; axis off;
	end
end
