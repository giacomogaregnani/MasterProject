function mesh = toy_micro_mesh_generator2( p, options )
%TOY_MICRO_MESH_GENERATOR1 Summary of this function goes here
%   Detailed explanation goes here

%% TODO: Make a very simple mesh generator! with very few DOF!
s= 0.5;
mesh = structured_mesh([-s,s,-s,s],3,struct('bc','periodic'));
mesh.node([5,6,8,9],:) = [ ...
	-s/2 - cos((p(2)+p(1))*2*pi/3)/8, -s/2 + cos((p(1)+p(2))*pi)/8;
	-s/2 - cos((p(2)+p(1))*2*pi/3)/8,  s/2 + sin((p(1)+p(2))*pi)/8;
	 s/2 + cos((p(2)+p(1))*2*pi/3)/8,  s/2 + cos((p(1)-p(2))*pi)/8;
	 s/2 + cos((p(2)+p(1))*2*pi/3)/8, -s/2 + sin((p(1)-p(2))*pi)/8];
mesh.elem(9:10,:) = [];
%mesh = label(mesh);
%mesh = bisect(mesh,'all');
end

