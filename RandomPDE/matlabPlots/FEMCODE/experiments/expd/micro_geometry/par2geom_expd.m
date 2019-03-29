function [bdmesh, bdgmsh] = par2geom_expd(par, varargin)
%PAR2GEOM_EXPB Summary of this function goes here
%   Detailed explanation goes here
holes = cell(1,1);
c=(1/3+1/8)/2;
s=(1/3-1/8)/2;
holes{1} = ...
	[-c+s*sin(2*pi*par+0*pi/4), -c+s*sin(-2*pi*par + 2*pi/4);
	  c+s*sin(2*pi*par+1*pi/4), -c+s*sin(2*pi*par + 3*pi/4);
	  c+s*sin(2*pi*par+4*pi/4),  c+s*sin(-2*pi*par + 6*pi/4);
	 -c+s*sin(2*pi*par+5*pi/4),  c+s*sin(2*pi*par + 7*pi/4)];

[bdmesh, bdgmsh]= geom2D_gen_polyg_holes(holes);
end

