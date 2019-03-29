function [bdmesh, bdgmsh] = micBdmesh_X(x, ~)
%MICBDMESH_X Summary of this function goes here
%   Detailed explanation goes here
par = x2par_X(x);

bdmesh.node = [...
  -1/2,0; 
  par(1),par(2); 
     0,-1/2; 
   1/2,-1/2; 
   1/2, 0; 
   1/2, 1/2; 
     0, 1/2; 
  -1/2, 1/2];
bdmesh.elem = [1,2; 2,3; 3,4; 4,5; 5,6; 6,7; 7,8; 8,1];


bdgmsh.node = bdmesh.node;
bdgmsh.line = bdmesh.elem;
bdgmsh.lineloop = {1:8};
bdgmsh.periodicline = [3, -6; 5, -8];
bdgmsh.planesurface = {1};
bdgmsh.box = [-0.5,0.5,-0.5,0.5];
end

