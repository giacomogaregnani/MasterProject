function mesh = toy_micro_mesh_generator4( p, options )
%TOY_MICRO_MESH_GENERATOR1 Summary of this function goes here
%   Detailed explanation goes here

%% TODO: Make a very simple mesh generator! with very few DOF!
s= 0.5;
mesh = structured_mesh([-s,s,-s,s],3);
mesh.node([6,10,11,7],:) = [ ...
	-s/2 - cos((p(2)+p(1))*2*pi/3)/8, -s/2 + cos((p(1)+p(2))*pi)/8;
	-s/2 - cos((p(2)+p(1))*2*pi/3)/8,  s/2 + sin((p(1)+p(2))*pi)/8;
	 s/2 + cos((p(2)+p(1))*2*pi/3)/8,  s/2 + cos((p(1)-p(2))*pi)/8;
	 s/2 + cos((p(2)+p(1))*2*pi/3)/8, -s/2 + sin((p(1)-p(2))*pi)/8;];
mesh.elem(9:10,:) = [];
%mesh = label(mesh);
%mesh = bisect(mesh,'all');


mesh.border{1,1} = 1;      % lower left
mesh.border{1,2} = 13;   % upper left
mesh.border{2,1} = 4; % lower right
mesh.border{2,2} = 16;   % upper right

mesh.border{3,1} = 1:4;   % down
mesh.border{3,2} = 13:16;        % up
mesh.border{1,3} = 1:4:13;             % left
mesh.border{2,3} = 4:4:16; % right

end

