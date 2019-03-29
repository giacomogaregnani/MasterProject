function [bdmesh, bdgmsh] = par2geom_expc(par, varargin)
%PAR2GEOM_EXPB Summary of this function goes here
%   Detailed explanation goes here

theta = par/3;
ang = 2*pi*par;
shifts = [0, 2*pi/3, 4*pi/3];
eyes = 0.145 + 0.035*sin(ang + shifts);

[bdmesh, bdgmsh] = geom2D_expc(eyes, theta);
end

