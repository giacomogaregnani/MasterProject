function bdmesh = cg2d_n1(x, epsilon)
%cg2d_n1 defines a porous medium in 2D
%   This is an m-file that describes 2d periodic porouse geometry. Given
%   macro coordinate x and the porosity size epsilon, it gives back the
%   geometry description of the fluid part in the area 
%   [x-epsilon,x+epsilon]^2 scaled to translated to the domain interval 
%   [-0.5, 0.5]^2.
%       
%       INPUT
%   x               coordinate vector [x1, x2] of a macro point
%   varargin        these functions usually use epsilon as second parameter
%
%       OUTPUT
%   bdmesh          is a polygon description.
%
%       GEOMETRY DESCRIPTION
%
%   f1-f8           8 points on the boundary
%
%   B--7--------------6--B       1/2
%   |X/                \X|
%   8                    5           X - solid part
%   |                    |
%   |      fluid         |
%   |      part          |            fluid parts
%   1                    4            are encircled
%   |X\                /X|            counterclockwise
%   B--2--------------3--B      -1/2
%
% -1/2                  1/2

if epsilon > 0.5
	error('epsilon is too small, definition blows up for epsilon -> 1');
end

% dependence of local coordinates in the scaled cell domain on the
% global coordinates
f{1} = @(y)[-1,-0.5+sin(y(1)+y(2))/5]/2;
f{2} = @(y)[-0.5+sin(y(1)-y(2))/5, -1]/2;
f{3} = @(y)[0.5+sin(y(1))/5, -1]/2;
f{4} = @(y)f{1}(y) + [1,0];
f{5} = @(y)[1, 0.5+sin(y(2))/5]/2;
f{6} = @(y)f{3}(y) + [0,1];
f{7} = @(y)f{2}(y) + [0,1];
f{8} = @(y)f{5}(y) - [1,0];
    
% find fixed points of the coordinates
if epsilon == 0
	numiter = 1;
else
	numiter = ceil(-20/log10(epsilon));
end
g=cell(8,1);
for i=1:8
	g{i} = [0,0];
	for j=1:numiter
		g{i} = f{i}(x + epsilon * g{i});		
	end
end

bdmesh.node = [g{1}; g{2}; g{3}; g{4}; g{5}; g{6}; g{7}; g{8}];

bdmesh.elem = [1,2; 2,3; 3,4; 4,5; 5,6; 6,7; 7,8; 8,1];

bdmesh.box=[-0.5, 0.5, -0.5, 0.5];

bdmesh.bdnode{1,1} = []; % lower left
bdmesh.bdnode{1,2} = []; % upper left
bdmesh.bdnode{2,1} = []; % lower right
bdmesh.bdnode{2,2} = []; % upper right

bdmesh.bdnode{3,1} = [2, 3]; % down
bdmesh.bdnode{3,2} = [7, 6]; % up
bdmesh.bdnode{1,3} = [1, 8]; % left
bdmesh.bdnode{2,3} = [4, 5]; % right

bdmesh.bdelem{3,1} = 2; % down
bdmesh.bdelem{3,2} = 6; % up
bdmesh.bdelem{1,3} = 8; % left
bdmesh.bdelem{2,3} = 4; % right
end

