function [bdmesh, bdgmsh] = cg2d_1(x, varargin)
%cg2d_1 gives a sample 2D cell geometry with one quadrilateral solid hole
%
%   This is an m-file that describes 2d periodic porouse geometry. Given
%   macro coordinate x it gives back the geometry description of the fluid 
%   part in the area [x-epsilon/2,x+epsilon/2]^2 scaled and translated 
%   to the interval [-0.5,0.5]^2.
%       
%       INPUT
%   x               coordinate vector [x1, x2] of a macro point
%   varargin        these functions usually use epsilon as second parameter
%
%       OUTPUT
%   bdmesh          is a polygon description.
%
%       GEOMETRY DESCRIPTION
%   The geometry has only one solid part that is strictly included in the
%   cell. The coordinates of the nodes in the output depend only on the 
%   center of the cell, that is x.
%
%  1/2-X------------------<-------X
%      |                          |
%     \|/       2---->-----3      |
%      |       /|\         |      |
%      | fluid  | solid    |      |
%      | part   |  part   \|/     |
%      |        |          |     /|\
%      |        1--__      |      |
%      |             -<-_  |      |
%      |   (clockwise    --4      |
%      |   direction)             |
% -1/2-X---------------->---------X
%      |                          |
%    -1/2                        1/2

holes = {[ ...
	-1/4 - cos((x(2)+x(1))*2*pi/3)/8, -1/4 + cos((x(1)+x(2))*pi)/8;
	-1/4 - cos((x(2)+x(1))*2*pi/3)/8,  1/4 + sin((x(1)+x(2))*pi)/8;
	 1/4 + cos((x(2)+x(1))*2*pi/3)/8,  1/4 + cos((x(1)-x(2))*pi)/8;
	 1/4 + cos((x(2)+x(1))*2*pi/3)/8, -1/4 + sin((x(1)-x(2))*pi)/8]};
 
 [bdmesh, bdgmsh]= geom2D_gen_polyg_holes(holes);
end

