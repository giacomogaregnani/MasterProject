function bdmesh = cg2d_n2(x, varargin)
%cg2d_n2 defines a porous medium in 2D
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
%   8          fluid     5           X - solid part
%   |       9    part    |
%   |      /X\           |
%   |    12 X 10         |            fluid parts
%   1      \X/           4            are encircled
%   |X\     11           /X|            counterclockwise
%   B--2--------------3--B      -1/2
%
% -1/2                  1/2


bdmesh.node = [...
	-0.5, -0.1; -0.1, -0.5;...
	 0.1, -0.5;  0.5, -0.1; ...
	 0.5,  0.1;  0.1,  0.5;   ...
	-0.1,  0.5; -0.5,  0.1;  ...
	 0,0.3; 0.3,0; 0, -0.3; -0.3,0];

bdmesh.elem = [1,2; 2,3; 3,4; 4,5; 5,6; 6,7; 7,8; 8,1;
	9,10; 10,11; 11,12; 12,9];

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