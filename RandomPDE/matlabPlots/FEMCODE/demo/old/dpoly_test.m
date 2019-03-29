%% 2D EXAMPLE
clear options meshf; options.inside = 1;
N = 500; 

pf=rand(N,2)*3-0.5; % random points
meshf.node = [0,0; 2,0; 2,1; 1,1; 1,2; 0,2];
meshf.elem = [1,2; 2,3; 3,4; 4,5; 5,6; 6,1];
[ ~, closestf, insidef, ~, ~] = dpoly(pf, meshf.node, meshf.elem, options);

figure; hold on; axis equal; % plotting
simpplot(meshf,struct('facecolor','green'));
for i=1:size(pf,1)   % plot projections
	if insidef(i), color = 'blue'; else color='red'; end
	line([pf(i,1), closestf(i,1)],[pf(i,2), closestf(i,2)], 'Color',color);
end
hold off;

%% 3D EXAMPLE

clear options mesh; options.inside = 1;
N=500;

mesh.node = [-0.5,-0.5,0; 1,0,0; 0,1,0; -1,-1,1; -1,-1,-1];
mesh.elem = [1,4,2; 2,4,3; 3,4,1; 2,5,1; 3,5,2; 1,5,3];
p3= rand(N,3)*3-1.5;
[~, closest3, inside3, ~, ~] = dpoly(p3, mesh.node, mesh.elem, options);

%% PLOTTING 3D
figure; hold on; axis equal;
simpplot(mesh,struct('color','yellow','facealpha',0.2));
for i=1:size(p3,1)   % plot projections
	if inside3(i), color = 'blue'; else color='red'; end
	line([p3(i,1), closest3(i,1)], [p3(i,2), closest3(i,2)], ...
		 [p3(i,3), closest3(i,3)], 'Color', color);
end
hold off;