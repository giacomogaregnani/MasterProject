function [bdmesh, bdgmsh] = geom2D_expc(hole_size, theta)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%% EYE 1
sc = 0.25;
hole_shift1 = sc*[-sqrt(3)/2,1/2];
hole_shift2 = sc*[sqrt(3)/2,1/2];
hole_shift3 = sc*[0,-1];
hole_shape =  [1,0; 1/2, -sqrt(3)/2; -1/2, -sqrt(3)/2; -1,0; ...
    -1/2, sqrt(3)/2; 1/2, sqrt(3)/2];
holes{1} = bsxfun(@plus,hole_shape*hole_size(1),hole_shift1);
holes{2} = bsxfun(@plus,hole_shape*hole_size(2),hole_shift2);
holes{3} = bsxfun(@plus,hole_shape*hole_size(3),hole_shift3);

ang = 2*pi*theta;
R = [cos(ang), -sin(ang); sin(ang), cos(ang)];

for i=1:numel(holes)
    holes{i} = holes{i}*R;
end

[bdmesh, bdgmsh]= geom2D_gen_polyg_holes(holes);
end

