function mesh = periodize(mesh, box)
%PERIODIZE Summary of this function goes here
%   Detailed explanation goes here

if isfield(mesh,'periodic') && mesh.periodic
	display('mesh already periodic'); 
	return; 
end

d = size(mesh.node,2);
dimensions = box(2:2:end) - box(1:2:end);
tol = 1e-8 * min(dimensions);
N = size(mesh.node,1);

%% Find boundary nodes
bdnode=zeros(N,1);
for i=1:d
	bdnode = bdnode + ...
		(abs(mesh.node(:,i) - box(2*i-1)) < tol) + ...
		(abs(mesh.node(:,i) - box(2*i)) < tol);
end
bdnode = find(bdnode);

%% Master-slave relationship
master = bdnode;
for i=1:numel(bdnode)
	x = abs( bsxfun(@minus, mesh.node(bdnode(1:i-1),:), mesh.node(bdnode(i),:)));
	x = x .* bsxfun(@minus, x, (box(2:2:end) - box(1:2:end)));
	x = find(sum(abs(x),2) < tol);
	if ~isempty(x)
		master(i) = master(x(end));
	end
end
idx = (master ~= bdnode);
master = master(idx);
slave = bdnode(idx);

%% SUBSTITUTION
%1. make substitution master -> slave
%2. shift node numbers to remove gaps
%3. remove unnecessary nodes
subs = 1:N;
subs(slave) = master;

idx = true(1,size(mesh.node,1));
idx(slave) = false;
subs(idx) = 1:sum(idx);
subs(~idx) = subs(subs(~idx));

mesh.elem = subs(mesh.elem);

mesh.periodic = true;
mesh.box = box;

mesh.node(slave,:) = [];
end

