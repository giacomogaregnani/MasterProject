function drawpolygon(id, poly, opt)
%DRAWTOFILE Summary of this function goes here
%   Detailed explanation goes here

if nargin < 3
  opt = '';
end

N = size(poly,1);
fprintf(id,'\\draw');
if ~isempty(opt)
  fprintf(id,'[%s]',opt);
end
  fprintf(id,' ');
for i=1:N
  fprintf(id,'(%f, %f) -- ', poly(i,1), poly(i,2));
end
fprintf(id,'cycle;\n');
end

