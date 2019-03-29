function [bdmesh, bdgmsh] = pormat_expe(x, varargin)
%CG2D_EXP Summary of this function goes here
%   Detailed explanation goes here

% THIS INFO IS ALSO IN ANOTHER SCRIPT
r = -x(:,1).^2/8 + (3-x(:,2))/3;

[bdmesh, bdgmsh] = par2geom_expe(r);
end

