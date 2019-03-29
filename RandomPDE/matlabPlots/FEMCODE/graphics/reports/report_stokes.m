function report_stokes(mesh, vp, ufemspace, usol, pfemspace, psol)
%REPORT_STOKES Summary of this function goes here
%   Detailed explanation goes here

NS = size(usol,3);
dim = size(mesh.elem,2)-1;
	
fprintf(1,'Velocity elements: %s \n',vp.elemtype);
fprintf(1,'Pressure elements: %s \n',vp.pelemtype);
fprintf(1,'DOF: %d \n',dim*ufemspace.ndof+pfemspace.ndof);
fprintf(1,'NT = Number of Elements: %d \n',size(mesh.elem,1));
fprintf(1,'N  = Number of Nodes: %d \n',size(mesh.node,1));

for i=1:NS
	fprintf(1,'Stokes problem solution, ID = %d.\n',i);

	%% Plot mesh
	clf;
	simpplot(mesh);
	title('mesh plot');
	drawnow; snapnow; hold off;

	%% Plot pressure
	clf;
	simpplot_sol(mesh,psol(:,1,i));
	title('pressure solution');
	drawnow; snapnow; hold off;
	
	%% plot velocity
	for j=1:dim
		clf;
		simpplot_sol(mesh,usol(:,j,i));
		title(['velocity solution (coor ' num2str(j) ')']);
		drawnow; snapnow; hold off;
	end
end
clf;
end

