function [mesh, idx, signed_volume] = fixorder(mesh)
%FIXORDER sets orientation of all elements counter-clockwise
% 
%   [mesh, idx, volume] = fixorder(mesh) computes signed measure of all
%   elements and switches the first two nodes in each elements where this 
%   measure is negative, i.e. it assures the orientation of all elements is 
%   counter-clockwise. Additional outputs are the indices of elements that
%   were changed (idx) and the signed volume (signed_volume).

signed_volume = simplex_volume(mesh, true);
idx = find(signed_volume < 0);
mesh.elem(idx, [1 2]) = mesh.elem(idx, [2 1]);
