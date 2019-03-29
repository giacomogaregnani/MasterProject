function filltofile( mesh, id, color, intervals )
%DRAWTOFILE Summary of this function goes here
%   Detailed explanation goes here
if nargin < 3 || isempty(color)
	color = '\colorb';
end
if size(mesh.elem,2) > 2
	for i=1:size(mesh.elem,1)
		fprintf(id,'\\filldraw[fill = %s] ',color);
		for j=1:size(mesh.elem,2)
			fprintf(id,'(%f, %f)', ...
				mesh.node(mesh.elem(i,j),1), ...
				mesh.node(mesh.elem(i,j),2));
				fprintf(id,' -- ');
		end
		fprintf(id,'cycle;\n');
	end
	return;
else
	edge = mesh.elem;
end
NE = size(edge,1);

if nargin<4, intervals = [1,NE]; end

for k=1:size(intervals,1)	
	fprintf(id,'\\filldraw[fill = %s] ',color);
	for i=intervals(k,1):intervals(k,2)
		fprintf(id,'(%f, %f) -- ', ...
			mesh.node(edge(i,1),1), mesh.node(edge(i,1),2));
	end
	fprintf(id,'cycle;\n');
end
end

