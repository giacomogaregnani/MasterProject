function bdmesh = get_micro_geometry( fg, sector, options )
%GET_MICRO_GEOMETRY3 Summary of this function goes here
%   Detailed explanation goes here
% TODO: update for noninclusive geometries

%% PARAMETER READING
if nargin<2, 
	error('sector must be specified');
end
if nargin<3, options = struct; end
if ~isfield(sector,'centre'), 
	error('sector centre not set, continuing with the origin');
end
dim = numel(sector.centre);
if ~isfield(fg,'epsilon')
	warning('WARNING:epsilonNotSet', ...
		'epsilon not set, continuint with 0.01');
	fg.epsilon = 0.01;
end
if ~isfield(sector,'delta')
	warning('WARNING:deltaNotSet', ...
		'delta not set, continuint with epsilon value');
	fg.delta = fg.epsilon;
end
if numel(fg.epsilon) == 1,
	fg.epsilon = repmat(fg.epsilon,[1,dim]);
end
if numel(sector.delta) == 1,
	sector.delta = repmat(sector.delta,[1,dim]);
end
if any(0.999*fg.epsilon > sector.delta)
	error('one must provide delta >= epsilon');
end
if ~isfield(options, 'base')
	options.base = 'origin';
end
if strcmp(options.base,'origin')
	shift = zeros(1,dim);
elseif strcmp(options.base,'shifted')
	shift = fg.epsilon/2;
end
if ~isfield(options,'periodic')
	options.periodic = false;
end

%% DETERMINING WHICH CELLS AFFECT OUR GEOMETRY
kmin = ceil((sector.centre+shift-sector.delta/2)./fg.epsilon - 1/2);
kmax =floor((sector.centre+shift+sector.delta/2)./fg.epsilon + 1/2);
geom=cell(kmax-kmin+1);
if (dim==3)	
	for k1=kmin(1) : kmax(1)
		for k2= kmin(2) : kmax(2)
			for k3 = kmin(3) : kmax(3)
				centre = shift + [k1,k2,k3].*fg.epsilon; 
				geom{k1-kmin(1)+1,k2-kmin(2)+1,k3-kmin(3)+1} = ...
					fg.handle(centre, fg.epsilon);  % local geometry
			end
		end
	end
elseif (dim == 2)
	for k1=kmin(1) : kmax(1)
		for k2= kmin(2) : kmax(2)
			centre = shift + [k1,k2].*fg.epsilon; 
			geom{k1-kmin(1)+1,k2-kmin(2)+1} = ...
				fg.handle(centre, fg.epsilon);  % local geometry
		end
	end
else
	error('Unsupported dimension');
end
clear centre;

NG = 0;
bdmesh.node = [];
bdmesh.elem = [];
gs = size(geom);
if (dim == 3)
	for i=1:gs(1)
	for j=1:gs(2)
	for k=1:gs(3)
		N = size(geom{i,j,k}.node, 1);
		% shift nodes
		geom{i,j,k}.node = geom{i,j,k}.node + repmat([i-1,j-1,k-1], [N,1]);
		geom{i,j,k}.offset = NG;
		% shift elements
		geom{i,j,k}.elem = geom{i,j,k}.elem + geom{i,j,k}.offset;
		% remove unnecessary boundary faces
		trm = [];
		if i>1, trm = [trm; geom{i,j,k}.bdelem{1,3,3}]; end
		if j>1, trm = [trm; geom{i,j,k}.bdelem{3,1,3}]; end
		if k>1, trm = [trm; geom{i,j,k}.bdelem{3,3,1}]; end
		if i<gs(1), trm = [trm; geom{i,j,k}.bdelem{2,3,3}]; end
		if j<gs(2), trm = [trm; geom{i,j,k}.bdelem{3,2,3}]; end
		if k<gs(3), trm = [trm; geom{i,j,k}.bdelem{3,3,2}]; end
		geom{i,j,k}.elem(trm,:)=[];
		bdmesh.node(NG+1:NG+N, 1:3) = geom{i,j,k}.node;
		bdmesh.elem =[bdmesh.elem; geom{i,j,k}.elem];
		NG = NG + N;
	end
	end
	end
elseif (dim == 2)
	for i=1:gs(1)
	for j=1:gs(2)
		N = size(geom{i,j}.node, 1);
		% shift nodes
		geom{i,j}.node = geom{i,j}.node + repmat([i-1,j-1],[N,1]);
		geom{i,j}.offset = NG;
		% shift elements
		geom{i,j}.elem = geom{i,j}.elem + geom{i,j}.offset;
		% remove unnecessary boundary edges
		trm = [];
		if i>1, trm = [trm; geom{i,j}.bdelem{1,3}]; end
		if j>1, trm = [trm; geom{i,j}.bdelem{3,1}]; end
		if i<gs(1), trm = [trm; geom{i,j}.bdelem{2,3}]; end
		if j<gs(2), trm = [trm; geom{i,j}.bdelem{3,2}]; end
		geom{i,j}.elem(trm,:)=[];
		bdmesh.node(NG+1:NG+N, 1:2) = geom{i,j}.node;
		bdmesh.elem =[bdmesh.elem; geom{i,j}.elem];
		NG = NG + N; 
	end
	end
else
	error('Unsupported dimension');
end

%% FIXING BOX AND PERIODIC
centre = (sector.centre - shift)./fg.epsilon - kmin;
bdmesh.box = [centre-sector.delta/fg.epsilon/2; 
	centre+sector.delta/fg.epsilon/2];
bdmesh.box = reshape(bdmesh.box,1,[]);
bdmesh.periodic = options.periodic;
end

