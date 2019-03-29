function [bdmesh, bdgmsh] = cg3d_1(x, varargin)
%cg3d_1 gives a sample 3D cell geometry with one tetrahedral solid hole
%          G---------H
%         /|        /|
%        / |       / |                Z  Y
%       E---------F  |                | /
%       |  C------|--D   AND 9 in     |/
%       | /       | /    the centre   O------X
%       |/        |/
%       A---------B		

tetvert = [...
	-0.4,-0.4,-0.4;
     0.4,-0.35, -0.4;
     0,0.4,-0.4;
     0.3 + sin(x(1))/10,  cos(x(2))/4, 0.15 + cos(x(3))/8];
[bdmesh, bdgmsh] = geom3D_gen_tet_hole( tetvert );
end

