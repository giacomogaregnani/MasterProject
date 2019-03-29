function mesh = remove_duplicate_nodes(mesh, options)
%REMOVE_DUPLICATE_NODES Summary of this function goes here
%   Detailed explanation goes here

if nargin < 2, options = struct; end
if ~isfield(options,'tol'), options.tol = 1e-5; end
if ~isfield(options,'N'),   options.N   = 10; end

d = size(mesh.node,2);

for i=1:options.N
  direc = rand(d,1);
  f = mesh.node * direc;
  [~, I] = sort(f);
  dists = sum((mesh.node(I(1:end-1),:)-mesh.node(I(2:end),:)).^2,2);
  cl = find(dists < options.tol^2);
  % substitute I(cl+1) to I(cl);
  subs = sort([I(cl), I(cl+1)],2);
  mesh.elem =subs_array(mesh.elem, subs(:,2), subs(:,1));
  mesh = renumber(mesh);
end

end

