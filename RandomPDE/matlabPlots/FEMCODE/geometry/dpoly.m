function [distance, closest, inside, whe, whn] = ...
        dpoly(points, node, elem, options)
%DPOLY Computes distances of a set of points to a boundary mesh
%
% [ distance, closest, inside, whe, whn ] = dpoly(points, node, elem) 
% returns:
%
%   distance = the squares of Euclidean distances (distance) of an 
%   NP x d - array of points in R^d to a dim-dimensional (dim <= d) mesh 
%   (node, elem). 
%
%   closest = NP x d array that gives coordinates of the closest points on 
%   bdmesh to which the distances are minimal.
%
%   inside = logical array of size NP x 1 specifying whether the points are
%   inside or outside. (true = inside)
%
%   whe = array of size NP x 1, the indices of the elements in bdmesh in 
%   which lie the closest points.
%
%   whn = logical array of size NP x (dim+1), where whn(i,:) is true at
%   those indices that define a minimal subsimplex in whe(i) to which the
%   distance from points(i) is the same as distance(i).
%
% [ ... ] = dpoly(points, node, elem, options) can specify additional
% options for the alogirthm.
%
%   options.inside = true / false : determines if we want to compute the
%   output inside (can save some time).

%% PARAMETER CHECKING
if (nargin < 3), error('Not enought parameters'); end
if (nargin < 4), options=struct; end
if ~isfield(options,'inside'),   options.inside = 0;   end % default values

NP   = size(points,1); % Number of Points
d    = size(points,2); % dimension of space (R^d)
dim  = size(elem,2)-1; % dimension of simplical mesh (bd)
BIG  = 10^20;          % big number, greater than any distance
tol  = 1e-7;           % tolerance for equality of sines of angles
tol2 = 1e-14;          % tolerance for equality of squared distances

%% INITIALIZE OUTPUT
distance = repmat(BIG, [NP,1]);
closest  = zeros(NP,d);
whe    = zeros(NP,1);
whn    = zeros(NP,dim + 1);
inside   = zeros(NP,1);

%% SIMPLE AND ERROR CASES
if (NP == 0), return; end
if (dim > d), error('Invalid dimensions'); end
if (dim + 1 < d) && options.inside, 
        error('No inside/outside defined by elem'); end
if (dim == 0) % node, elem are only points
	for i=1:size(elem,1)
		new_distance = sum((difference(points,node(elem(i),:))).^2,2);
		idx = new_distance < distance;
		distance(idx) = new_distance(idx);
		whe(idx) = i;
	end
        closest = node(elem(whe),:);
        whn   = ones(NP,1);
		return;
end

%% GENERAL CASE
new_distance = zeros(NP,1);         % distance to elem(i,:)
new_closest  = zeros(size(points)); % closest point to elem(i,:)
new_wnode    = zeros(NP,dim+1);     % wnode field to elem(i,:)
for i=1:size(elem,1)                % process one elem at a time
    % project points on the subspace spanned by node(elem(i,:),:)
    [proj, inface] = project(points,node(elem(i,:),:)); 
    new_distance(inface) = ...
        sum(difference(points(inface,:), proj(inface,:)).^2, 2);
    new_closest(inface,:)= proj(inface,:);
    new_wnode(inface,:)  = 1;
     
    % recursively call this function with all the faces of elem(i,:)
    [~, new_closest(~inface,:), ~, sub_welem, sub_wnode] = ...
        dpoly(points(~inface,:), node, ...
        get_subsimplices(elem(i,:),'opposite'));
    new_distance(~inface) = ...
        sum(difference(points(~inface,:), new_closest(~inface,:)).^2,2);
    new_wnode(~inface,:) = 0; 
    for k=1:dim+1
        new_wnode(logjoin(~inface, sub_welem == k), [1:k-1, k+1:dim+1])=...
			sub_wnode(sub_welem == k,:);
	end
        
	% difference between distance to elem(i,:) and minimal to elem(1:i-1,:)
	subt = sum(distance.^2,2) - sum(new_distance.^2,2);
        
	if ~options.inside % simple case - no inside / outside computation
		update = subt > 0;
	else
		update = subt >= tol2;     % move to update
		ver = abs(subt) < tol2;    % more verification
                
		% 1. greater dimension of subsimplex is preffered
		aux = sum(new_wnode(ver,:), 2) - sum(whn(ver,:), 2);
		update(ver) = aux > 0;  % move to update
		ver(ver) = aux == 0;    % more verification
                
        % 2. subsimplices of equal dimension have complete preference order
        NA   = sum(ver);
		ibd1 = sort(whn(ver,:) .* elem(whe(ver),:),2);
		ibd2 = sort(new_wnode(ver,:) .* repmat(elem(i,:),[NA,1])  ,2);
		subt = sum((ibd1-ibd2) .* repmat((NP+1).^(dim:-1:0), [NA,1]), 2);
		update(ver) = subt >  0;  % move to update
		ver(ver)    = subt == 0;  % further verification
                
		% 3. for equal subsimplices, angle considerations must be done
		ang1 = get_angle(node, elem, points(ver,:), closest(ver,:), ...
			whe(ver,:), whn(ver,:));
		ang2 = get_angle(node, elem, points(ver,:), closest(ver,:), ...
			i, new_wnode(ver,:));
		update(ver) = maj(ang1-ang2, tol); % move to update
	end
	
	% update output
	distance(update) = new_distance(update);
	closest(update,:) = new_closest(update,:);
	whe(update) = i;
	whn(update,:) = new_wnode(update,:);
end

%% Detect inside / outside points
if options.inside
	mat = zeros(NP, dim+1, dim+1); % NP matrices dim times dim
	mat(:,:,dim + 1) = points - closest;
	for k=1 : dim
		mat(:,:,k) = node(elem(whe, k+1), :) - node(elem(whe, k), :);
	end
	inside = detn(mat) > 0;  % negative orientation => inside
end
end

%% ANGLES COMPUTATIONS
function output = get_angle(node, elem, apoints, aclosest, awelem, awnode)
if numel(awelem) == 1, awelem = repmat(awelem,[size(apoints,1),1]); end
N   = size(apoints,1);
NT  = size(elem,1);
dim = size(elem,2) - 1;
d   = size(node,2);
BIG = 10^20;

output = repmat(BIG, [N, dim+1]);
ip = ones(N,1);
for s = 1:dim
	id = (sum(awnode,2) == s);
	Nk = sum(id);
	if (Nk == 0), continue; end
	ld = decode(awnode(id,:), s);
	lp = decode(1-awnode(id,:), dim+1-s);
	vecd = zeros(Nk,d,s-1);
	vecp = zeros(Nk,d,dim+1-s);
	for r=1:s-1
		vecd(:,:,r) = node(elem( (ld(:,r+1)-1) * NT + awelem(id) ),:) ...
			- node(elem( (ld(:,1)  -1) * NT + awelem(id) ),:);
	end
	for m=1:dim+1-s
		vecp(:,:,m) = node(elem( (lp(:,m)-1) * NT + awelem(id) ),:) ...
			- node(elem( (ld(:,1)-1) * NT + awelem(id) ),:);
	end
	for r=1:s-1
		for m=1:dim+1-s
			alpha =  - dot(vecd(:,:,r), vecp(:,:,m), 2) ...
				./ dot(vecd(:,:,r), vecd(:,:,r), 2);
			vecp(:,:,m) = vecp(:,:,m) + repmat(alpha,[1,d]) .* vecd(:,:,r);
		end
	end
	out2 = zeros(Nk, dim+1-s);
	for m=1:dim+1-s
		vecp(:,:,m) = vecp(:,:,m) ./ ...
			repmat(sqrt(sum(vecp(:,:,m).^2,2)), [1,d]);
		out2(:,m) = abs(dot(apoints(id,:)-aclosest(id,:), vecp(:,:,m), 2));
	end
	[out2, ind2] = min(out2,[],2);
	ind3 = lp((1:Nk)'+Nk*(ind2-1));
	ind4 = id .* (1:N)';
	ind4 = ind4(ind4>0);
	awnode(ind4+N*(ind3-1)) = 1;
	output(ind4 + N*(ip(id)-1)) = out2;  % update output
	ip(id) = ip(id) + 1;                % update ip
end
end

%% find majorantly positive vectors
function ind = maj(mat, tol)
ind  = true(size(mat,1),1);
ind2 = ind;
for i=1:size(mat,2)
	ind(ind2)  = mat(ind2,i) > tol;
	ind2(ind2) = abs(mat(ind2,i)) <= tol;
end
end

%% Logical join
function log1 = logjoin(log1,log2)
    log1(log1) = log2;
end

%% find k-th 1 in an array of vectors
function list = decode(mat, k)
M= size(mat,1);
N= size(mat,2);
list = zeros(M,k);
for l=1:k
	list(:,l) = N - max(mat .* repmat(N-1:-1:0, [M,1]),[],2);
	mat((1:M) +M*(list(:,l)-1)') = 0;
end
end

%% Difference of two vectors (may be arrays of vectors)
function output=difference(vec1, vec2)
if numel(vec1) < numel(vec2)
	output = repmat(vec1, [size(vec2,1),1]) - vec2;
elseif numel(vec1) > numel(vec2)
	output = -repmat(vec2, [size(vec1,1),1]) + vec1;
else
	output = vec1-vec2;
end
end