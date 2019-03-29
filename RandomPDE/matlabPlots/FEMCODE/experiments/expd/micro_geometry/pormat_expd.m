function [bdmesh, bdgmsh] = pormat_expd(x, varargin)
%CG2D_EXP Summary of this function goes here
%   Detailed explanation goes here

r = x2par_expd(x);

[bdmesh, bdgmsh] = par2geom_expd(r);
end

