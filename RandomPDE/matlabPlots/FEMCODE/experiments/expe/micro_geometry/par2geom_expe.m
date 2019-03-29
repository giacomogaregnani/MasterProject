function [bdmesh, bdgmsh] = par2geom_expe(par, varargin)
%PAR2GEOM_EXPB Summary of this function goes here
%   Detailed explanation goes here
holes = cell(1,1);

% Rectangular hole, rotated by pi*par
R = [cos(pi*par), -sin(pi*par); sin(pi*par), cos(pi*par)]; 
holes{1} = (R * [0.3,0.15; 0.3,-0.15; -0.3,-0.15; -0.3,0.15]')';

[bdmesh, bdgmsh]= geom2D_gen_polyg_holes(holes);
end

