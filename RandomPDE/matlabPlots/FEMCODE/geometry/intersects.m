function onbd  = intersects(mesh, whe, bdmesh, options)
%INTERSECTS find elements of a mesh who have intersection with bdmesh
%
% onbd  = intersects( mesh, 'all', bdmesh) finds all elements of mesh that
% have an intersection with bdmesh and represents the output as a logical
% array.
%
% onbd  = intersects( mesh, whe, bdmesh) does the same but it works only
% with the elements with indices in whe.
%
% onbd  = intersects(mesh, whe, bdmesh, options) can be used to specify
% some advanced options:
%
%    options.approx = true / false : whether or not to use a faster bu only
%    approximate algorithm.


%% INITIALIZATION
if ~strcmp(whe,'all'), mesh.elem = mesh.elem(whe,:); end
if nargin < 4, options = struct; end
if ~isfield(options,'approx'), options.approx = false; end
d  = size(mesh.node,2);
d1 = size(mesh.elem,2) - 1;
d2 = size(bdmesh.elem,2) - 1;

NT = size(mesh.elem, 1);
NB = size(bdmesh.elem,1);

mesh.bary = get_rc(mesh);

%% FAST APPROXIMATE VERSION
if options.approx
	rad = sum((get_rc(mesh, [], 'all', 1, [], ...
		struct('adjust_to','bary')) - mesh.bary).^2,2);
	for j=2:d1+1
		rad = sum((get_rc(mesh,[], 'all', j, [], ...
			struct('adjust_to','bary')) - mesh.bary).^2,2);
	end
	distance = dpoly( mesh.bary, bdmesh.node, bdmesh.elem );
	onbd = distance < rad;
	return;
end

%% PRECISE COMPUTATION INITIALIZATION
comb = combnk(1:(d1+1)+(d2+1),d1+d2-d);
active = true([NT,1]);

%% COUPLE (A, rhs) INITIALIZATION
% TODO: IF NT*NB is too big, make a cycle!
A = zeros(NT*NB,d,d1+d2);
for s=1:d1
	A(:,:,s) = repmat(- get_vec(mesh,[],'all',[1,s+1]), [NB,1]);
end
for s=1:d2
	A(:,:,s+d1) = repval(bdmesh.node(bdmesh.elem(:,s+1),:) - ...
		bdmesh.node(bdmesh.elem(:,1),:),[NT,1]);
end
rhs = repval(-bdmesh.node(bdmesh.elem(:,1),:),[NT,1]) + ...
	repmat(get_rc(mesh, [], 'all', 1, [], ...
	struct('adjust_to','bary')), [NB,1]);

%% CYCLE THROUGH POSSIBLE "BOUNDARY CONDITIONS"
for k=1:size(comb,1)
	is1 = sum(comb(k,:) == 1) == 1;
	is2 = sum(comb(k,:) == 2) == 1; 
	free = setdiff(1:d1+d2,comb(k,comb(k,:)>2)-2);
	if is1,
		f1 = setdiff(1:d1,comb(k,:) - 2);
		f1 = f1(1);
		free = setdiff(free,f1);
	end
	if is2,
		f2 = setdiff(1:d2,comb(k,:) - 2 - d1);
		f2 = f2(1) + d1;
		free = setdiff(free,f2);
	end
	CA = A;
	Crhs = rhs; 
	if is1,
		for s=[1:f1-1, f1+1:d1+d2]
			CA(:,:,s)=CA(:,:,s) - CA(:,:,f1);
		end
		Crhs = rhs - CA(:,:,f1);
	end
	if is2,
		for s=[1:f2-1, f2+1:d1+d2]
			CA(:,:,s)=CA(:,:,s) - CA(:,:,f2);
		end
		Crhs = rhs - CA(:,:,f2);
	end		
	CA = CA(:,:,free); 
	sol = zeros(NT*NB, d1+d2);
	sol(:,free) = solve_all(CA, Crhs);
	if is1
		sol(:,f1) = 1 - sum(sol(:,1:d1),2);
	end
	if is2
		sol(:,f2) = 1 - sum(sol(:,d1+1:d2),2);
	end
	active2 = ((sum( sol(:, 1 : d1) < 0, 2) > 0) +  ...
		       (sum( sol(:, 1 : d1),     2) > 1) + ...
			   (sum( sol(:, d1+1 : d1+d2) < 0, 2) > 0) + ...
			   (sum( sol(:, d1+1 : d1+d2)    , 2) > 1)) > 0;
	active2 = sum(reshape(active2,[NT,NB]),2) == NB;
	active(active) = active2;
	A = A(repmat(active2,[NB,1]),:,:);
	rhs = rhs(repmat(active2,[NB,1]),:,:);
	mesh.elem = mesh.elem(active2,:);
	mesh.bary = mesh.bary(active2,:);
	NT = size(mesh.elem,1);
end

onbd = ~active;
end