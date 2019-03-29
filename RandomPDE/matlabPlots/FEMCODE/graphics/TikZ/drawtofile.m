function drawtofile( mesh, id, opt )
%DRAWTOFILE Summary of this function goes here
%   Detailed explanation goes here
if nargin < 3, opt = ''; end
if ~isempty(opt),
	opt = ['[' opt ']'];
end

if isfield(mesh,'periodic') && mesh.periodic
	mesh = deperiodize(mesh);
end
if size(mesh.elem,2) > 2
	edge = auxstructure(mesh,'sim1');
else
	edge = mesh.elem;
end
NE = size(edge,1);

for i=1:NE
	fprintf(id,'\\draw%s (%f, %f) -- (%f, %f);\n', opt,...
		mesh.node(edge(i,1),1), mesh.node(edge(i,1),2), ...
		mesh.node(edge(i,2),1), mesh.node(edge(i,2),2));
end
end

