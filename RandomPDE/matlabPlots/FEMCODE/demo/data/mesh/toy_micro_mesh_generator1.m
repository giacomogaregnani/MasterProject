function mesh = toy_micro_mesh_generator1(p, options)
%TOY_MICRO_MESH_GENERATOR1 Summary of this function goes here
%   Detailed explanation goes here

if nargin<2, options = struct; end
if ~isfield(options,'h'), options.h =0.1; end
if ~isfield(options,'maxit'), options.maxit = 45; end

%SPECIAL PARAMETERS FOR CDOMAIN
s = 0.5;
bdmesh.box = [-s,s,-s,s];
bdmesh.node = [-s,-s; -s,s; s,s; s,-s; 
	-s/2 - cos((p(2)+p(1))*2*pi/3)/8, -s/2 + cos((p(1)+p(2))*pi)/8;
	-s/2 - cos((p(2)+p(1))*2*pi/3)/8,  s/2 + sin((p(1)+p(2))*pi)/8;
	 s/2 + cos((p(2)+p(1))*2*pi/3)/8,  s/2 + cos((p(1)-p(2))*pi)/8;
	 s/2 + cos((p(2)+p(1))*2*pi/3)/8, -s/2 + sin((p(1)-p(2))*pi)/8;];
bdmesh.elem = [1,2; 2,3; 3,4; 4,1; 5,8; 8,7; 7,6; 6,5];
h = options.h;
	
is_wrong = true;
while is_wrong
	ns = round(2*s/h);
	pfix = [
		repmat(-s,[ns,1]), (-s:h:s-h)';
		(-s:h:s-h)', repmat(s,[ns,1]);
		repmat(s,[ns,1]), (s:-h:-s+h)';
		(s:-h:-s+h)', repmat(-s,[ns,1]) 	 
		];
	for i=5:size(bdmesh.elem,1)
		% force nodes between bdmesh.node(i,:) and bdmesh.node(i+1,:)
		n1 = bdmesh.node(bdmesh.elem(i,1),:); 
		n2 = bdmesh.node(bdmesh.elem(i,2),:); 
		dist = sqrt(sum((n2-n1).^2));
		numstep = ceil(dist/h);
		step = (n2-n1)/numstep;
		pfix = [pfix; repmat(n1,[numstep,1]) + ...
			repmat((0:numstep-1)',[1,2]) .* repmat(step, [numstep,1])];
	end
	pstart = [repmat((-s+h:h:s-h)', [ns-1,1]), repval((-s+h:h:s-h)',[ns-1,1])];
	mesh = mesh_polyhedron(bdmesh, options, pfix, pstart);
	N = size(mesh.node,1);
	subs =[1:ns, 1, ns+1:2*ns-1, 1, ns:-1:1, 2*ns-1:-1:ns+1, 2*ns: N-(2*ns+1)];
	subs2 = [ns+1, 2*ns+1:4*ns];
	mesh.elem = subs(mesh.elem);
	mesh.node(subs2,:)=[];
	mesh.box = [-s,s,-s,s];
	mesh.periodic = true;
	mesh = label(mesh);
	mesh = fixorder(mesh);
	
	%% TEST IF IS RIGHT
	qual = simplex_quality(mesh);
	allelem = unique(mesh.elem(:));
	if all(qual>0.1) && (numel(allelem) == size(mesh.node,1)) && ...
			(max(allelem) == size(mesh.node,1))
		is_wrong = false;
	end
	if is_wrong
		fprintf('Mesh generation failed, regenerating once more\n');
		h=h/sqrt(2);
		options.h = h;
		options.maxit = options.maxit + 10;		
	end
end

