function [fsol, ffemspace] = interpolate(cmesh, cfemspace, csol, fmesh, father, ffemspace)
%INTERPOLATE Interpolates a coarse FE function to a finer FE space
if nargin <= 5
  ffemspace = get_femspace(fmesh, cfemspace.elemtype);
end
s = size(csol);
fsol = zeros([ffemspace.ndof, s(2:end)]);
for i=1:ffemspace.ldof
  loc_lam = get_lc(cmesh, fmesh, father, ffemspace.nodelambda(i,:));
  fsol(get_dof(fmesh,'all',ffemspace,i),:,:,:) = ...
    evalf(cmesh, father, cfemspace, loc_lam, csol);
end
end