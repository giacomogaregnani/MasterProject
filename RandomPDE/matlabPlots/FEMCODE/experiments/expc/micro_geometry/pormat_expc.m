function [bdmesh, bdgmsh] = pormat_expc(x, varargin)
%CG2D_EXP Summary of this function goes here
%   Detailed explanation goes here

r = x2par_expc(x);

[bdmesh, bdgmsh] = par2geom_expc(r);
end

