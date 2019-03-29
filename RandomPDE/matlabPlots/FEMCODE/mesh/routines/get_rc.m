function rc = get_rc(mesh, sub, whe, whi, lambda, options)
%GET_RC returns RC (real coordinates) of points described using baryc. coor
% 
% rc = get_rc(mesh) returns barycenters of all the elements in a mesh.
%
% rc = get_rc(mesh, sub) returns barycenters of all simplices described by
% node indices in rows of the matrix sub. This may give nonsense results
% for periodic meshes if not used without understanding metrics on a torus.
% WE SET: [N1, N2] = size(sub);
% 
%   rc(i,:) contains the barycenter of nodes with indices sub(i,:). 
%
% rc = get_rc(mesh, [], whe, whi) returns barycenters of all simplices
% for each element in whe, where the vertices are described as local
% indices in the array whi.
% WE SET: N1 = numel(whe), N2 = numel(whi).
%
%   if (size(whi,1) == 1): rc(i,:) is the barycenter of nodes
% mesh.elem(whe(i), whi(:)). It is equivalent to: 
% get_rc(mesh, mesh.elem(whe, whi))
%
%   if (size(whi,1) >  1): rc(i,:) is the barycenter of nodes
% mesh.elem(whe(i), whi(i,:))
%
% CONVENTION: whe = [] or whe = 'all' implies whe = 1:NT
%             whi = [] implies whi = 1:dim+1
% 
% rc = get_rc(mesh, sub, whe, whi, lambda) takes the weights lambda and 
% computes the weighted barycenters. If (size(lambda,1) == 1), then we use 
% weights lambda(:) for every simplex. If (size(lambda,1) == N1), then in 
% computation of rc(i,:) we use the weights lambda(i,:).
%
% rc = get_rc(mesh, sub, whe, whi, lambda, options) can specify in the last
% parameter the method that we use to finally align the real coordinates. 
% This is available only if the mesh is periodic.
%
%   if opt.adjust_to == 'global' (default): rc(i,:) is aligned to lie in the
% box defined by mesh.box (if the mesh is periodic).
%
%   if opt.adjust_to = 'bary': rc(i,:) is adjusted to mesh.bary(whe(i),:)
% which is adjusted to lie in the box defined by the mesh.box

if ~isfield(mesh,'periodic'), mesh.periodic = false; end
if nargin < 2, sub = []; whe = 'all';
elseif nargin == 2, whe = [];
end
if nargin < 4, whi = []; end
if nargin < 5, lambda = []; end
if (nargin < 6) || ~isfield(options,'adjust_to')
	options.adjust_to = 'global'; % or 'bary'
end

NT = size(mesh.elem,1);
dim = size(mesh.elem,2) - 1;
d = size(mesh.node,2);

% INPUT PARSING
if isempty(sub)
	if strcmp(whe,'all'), 
		N1 = NT; 
	else
		N1 = size(whe,1); 
	end
	if isempty(whi), 
		whi = 1:dim+1; 
	end
	N2 = size(whi,2);
else
	[N1, N2] = size(sub);
end
if isempty(lambda), 
	lambda = repmat(1/N2, [1, N2]); 
end

% PROESSING
rc=zeros(N1,d);
if ~isempty(sub)
	for i=2:N2
		if size(lambda,1) == 1
			rc = rc + lambda(i) * get_vec(mesh, sub(:,[1, i]));
		else
			rc = rc + repmat(lambda(:,i),[1,d]) .* ...
				get_vec(mesh, sub(:,[1, i]));
		end
	end
	rc = rc + mesh.node(sub(:,1),:);
else
	for i=2:N2
		if size(lambda,1) == 1
			rc = rc + lambda(i) * get_vec(mesh, [], whe, whi(:,[1,i]));
		else
			rc = rc + repmat(lambda(:,i),[1,d]) .* ...
				get_vec(mesh, [], whe, whi(:,[1,i]));
		end
	end
	if size(whi,1) == 1
		if strcmp(whe,'all')
			rc = rc + mesh.node(mesh.elem(:,whi(1)),:);
		else
			rc = rc + mesh.node(mesh.elem(whe,whi(1)),:);
		end
	else
		if strcmp(whe,'all')
			rc = rc + mesh.node(mesh.elem(  (1:NT)'+NT*(whi(:,1)-1)  ),:);
		else
			rc = rc + mesh.node(mesh.elem(  whe+NT*(whi(:,1)-1)  ),:);
		end
	end
end

if mesh.periodic
	if strcmp(options.adjust_to,'global')
		rc=per_crop(rc,mesh.box);
	elseif strcmp(options.adjust_to,'bary')
		if ~isfield(mesh,'bary')
			bary = get_rc(mesh, [], whe);
			rc = bary + per_diff(rc-bary, mesh.box);
		else
			if strcmp(whe,'all')
				rc = mesh.bary + per_diff(rc-mesh.bary, mesh.box);
			else
				rc = mesh.bary(whe,:) + ...
					per_diff(rc-mesh.bary(whe,:), mesh.box);
			end
		end
	else
		error('UNKNOWN adjust_to OPTION');
	end
end

