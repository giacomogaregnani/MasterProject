function [bdmesh, bdgmsh] = cg3d_3(x, epsilon)
%cg3d_3 gives a sample 3D cell geometry with a connected hole
%          G---------H
%         /|        /|
%        / |       / |                Z  Y
%       E---------F  |                | /
%       |  C------|--D                |/
%       | /       | /                 O------X
%       |/        |/
%       A---------B	

%% INITIALIZATION
if nargin<2
    epsilon = 0;
end

%% NODES
v1 = @(xx)([0.2*sin(pi*xx(1)) + 0.3 ,0,0]);
v2 = @(xx)([0,0.2*sin(pi*xx(1)) + 0.3 ,0]);
v3 = @(xx)([0,0,0.2*sin(pi*xx(1)) + 0.3]);

[bdmesh, bdgmsh]= geom3D_connected2(v1, v2, v3, x, epsilon);


end

