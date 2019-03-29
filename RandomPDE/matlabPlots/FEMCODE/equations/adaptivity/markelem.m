function markedElem = markelem(res, options)
% MARK mark element.
%
% marks a subset of elements by Dorfler
% marking strategy. It returns an array of indices of marked elements
% markedElem such that sum(eta(markedElem)^2) > theta*sum(eta^2).

if nargin < 2, options = struct; end
if ~isfield(options,'method'), options.method = 'L2'; end
if ~isfield(options,'theta'),  options.theta = 0.5; end
if ~isfield(options,'unify'), options.unify = true; end

NT = size(res,1);

switch upper(options.method)
  case 'MAX'
    markedElem = bsxfun(@ge, res, options.theta * max(res, [], 1));
  case 'L2'
    NF = size(res,3);
    markedElem = false(NT,1,NF);
    for i=1:NF
      [sortedEta, idx] = sort(res(:,:,i), 1, 'descend');
      x = cumsum(sortedEta, 1);
      markedElem(idx(x < options.theta*x(end)),1,i) = true;
      markedElem(idx(1),1,i) = true; % at least one marked
    end
end

if options.unify
  markedElem = sum(markedElem,3) > 0;
end
