function sample = generate_sample(param, options)
%GENERATE_SAMPLE Summary of this function goes here
%   param.min
%   param.max
%   param.test
%
%   options.Npoints
%   options.regular
if nargin < 2, options = struct; end
if ~isfield(options,'Npoints'),   options.Npoints = 1000; end
if ~isfield(options,'regular'),   options.regular = false; end

P = numel(param.min);
if options.regular
  if numel(options.Npoints) == 1
    npts = repmat(round(options.Npoints^(1/P)), [1,P]);
  else
    npts = options.Npoints;
  end
  sample = zeros(1,0);
  for p=1:P
    sample = [repmat(sample,[npts(p),1]), ...
      repval(linspace(param.min(p),param.max(p),npts(p))', [size(sample,1),1])];
  end
else
  sample = bsxfun(@plus, ...
    bsxfun(@times,rand(prod(options.Npoints(:)),P), param.max-param.min), param.min);
end
sample = sample(param.test(sample),:); % leave only those in the domain
end

