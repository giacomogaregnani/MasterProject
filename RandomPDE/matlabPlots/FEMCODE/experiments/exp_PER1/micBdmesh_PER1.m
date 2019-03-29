function [bdmesh, bdgmsh] = micBdmesh_PER1(x, ~)
%MICBDMESH_X Summary of this function goes here
%   Detailed explanation goes here
par = x2par_X(x);

bdmesh.node = [...
  -1/2,-1/2; 
   1/2,-1/2; 
   1/2, 1/2; 
  -1/2, 1/2];

v = linspace(0,2*pi,20)';
v(end) = 0;
bdmesh.node = [bdmesh.node; [v, [v(2:end); v(1)]]];
N = numel(v);
mylist = 4+(1:N)';
bdmesh.elem = [1,2; 2,3; 3,4; 4,5; [mylist, [mylist(2:end); mylist(1)]]];


bdgmsh.node = bdmesh.node;
bdgmsh.line = bdmesh.elem;
bdgmsh.lineloop = {[1:8]};
bdgmsh.periodicline = [3, -6; 5, -8];
bdgmsh.planesurface = {1};
bdgmsh.box = [-0.5,0.5,-0.5,0.5];
end

