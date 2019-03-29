function drawatpoints(id, points, command )
%DRAWTOFILE Summary of this function goes here
%   Detailed explanation goes here

N = size(points,1);
for i=1:N
  fprintf(id,'\\begin{scope}[shift = {(%f, %f)}] %s \\end{scope}\n', points(i,1), points(i,2), command);
end
end

