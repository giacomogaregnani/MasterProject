function [bdmesh, bdgmsh]= geom3D_connected1( v1, v2, v3, x, epsilon )
%GEOM_3D_CONNECTED1 Summary of this function goes here
%   Detailed explanation goes here

bdmesh = [];

V =  [-0.5,-0.5,-0.5;
 0.5,-0.5,-0.5;
-0.5, 0.5,-0.5;
 0.5, 0.5,-0.5;
-0.5,-0.5, 0.5;
 0.5,-0.5, 0.5;
-0.5, 0.5, 0.5;
 0.5, 0.5, 0.5];

VI = [...
 -v1(x)-v2(x)-v3(x);
  v1(x)-v2(x)-v3(x);
 -v1(x)+v2(x)-v3(x);
  v1(x)+v2(x)-v3(x);
 -v1(x)-v2(x)+v3(x);
  v1(x)-v2(x)+v3(x);
 -v1(x)+v2(x)+v3(x);
  v1(x)+v2(x)+v3(x)];

x1 = x + epsilon * [-0.5, 0, 0];
x2 = x + epsilon * [ 0.5, 0, 0];

VII(:,:,1) = [...
  -v2(x1)-v3(x1);
  -v2(x2)-v3(x2);
  v2(x1)-v3(x1);
  v2(x2)-v3(x2);
  -v2(x1)+v3(x1);
  -v2(x2)+v3(x2);
  v2(x1)+v3(x1);
  v2(x2)+v3(x2)];

x1 = x + epsilon * [0, -0.5, 0];
x2 = x + epsilon * [0,  0.5, 0];

VII(:,:,2) = [...
 -v1(x1)-v3(x1);
  v1(x1)-v3(x1);
 -v1(x2)-v3(x2);
  v1(x2)-v3(x2);
 -v1(x1)+v3(x1);
  v1(x1)+v3(x1);
 -v1(x2)+v3(x2);
  v1(x2)+v3(x2)];

x1 = x + epsilon * [0, 0, -0.5];
x2 = x + epsilon * [0, 0,  0.5];

VII(:,:,3) = [...
 -v1(x1)-v2(x1);
  v1(x1)-v2(x1);
 -v1(x1)+v2(x1);
  v1(x1)+v2(x1);
 -v1(x2)-v2(x2);
  v1(x2)-v2(x2);
 -v1(x2)+v2(x2);
  v1(x2)+v2(x2)];

for j=1:3
    VII(:,j,j) = V(:,j);
end

bdgmsh.node = [V; VI; VII(:,:,1); VII(:,:,2); VII(:,:,3)];

%% LINES
bdgmsh.line = [1,2; 3,4; 5,6; 7,8; 
	           1,3; 2,4; 5,7; 6,8;
			   1,5; 2,6; 3,7; 4,8];
bdgmsh.line = [bdgmsh.line; bdgmsh.line + 8];
bdgmsh.line = [bdgmsh.line;
    17,9; 10,18; 19,11; 12,20; %% first group
    21,13; 14,22; 23,15; 16,24;
    17,19; 18,20; 21,23; 22,24;
    17,21; 18,22; 19,23; 20,24;
    25,26; 27,28; 29,30; 31,32; %% second group
    25,9; 11,27; 26,10; 12,28;
    29,13; 15,31; 30,14; 16,32;
    25,29; 26,30; 27,31; 28,32;
    33,34; 35,36; 37,38; 39,40; %% third group
    33,35; 34,36; 37,39; 38,40;
    33,9; 13,37; 34,10; 14,38;
    35,11; 15,39; 36,12; 16,40];
bdgmsh.lineloop = {[5, 11, -7, -9]; % ACGE
	               [6, 12, -8, -10]; % BDHF	 
	               [1, 10, -3, -9]; % ABFE
	               [2, 12, -4, -11]; % CDHG
                   [1, 6, -2, -5]; %ABDC
	               [3, 8, -4, -7]; %EFGH
                   [33,39,-35,-37]; %% HOLES %%
                   [34,40,-36,-38];
                   [41,54,-43,-53];
                   [42,56,-44,-55];
                   [57,62,-58,-61];
                   [59,64,-60,-63];
				   [25,17,-27,-33]; %% 1
                   [29,19,-31,-35];
                   [25,21,-29,-37];
                   [27,23,-31,-39];
                   [26,38,-30,-22]; %% 2
                   [28,40,-32,-24];
                   [26,34,-28,-18];
                   [30,36,-32,-20];
                   [45,21,-49,-53]; %% 3
                   [47,22,-51,-54];
                   [41,47,-13,-45];
                   [43,51,-15,-49];
                   [46,55,-50,-23]; %% 4
                   [48,56,-52,-24];
                   [14,48,-42,-46];
                   [16,52,-44,-50];
                   [61,69,-17,-65]; %% 5
                   [62,71,-18,-67];
                   [57,67,-13,-65];
                   [58,71,-14,-69];
                   [19,70,-63,-66]; %% 6
                   [20,72,-64,-68];
                   [15,68,-59,-66];
                   [16,72,-60,-70]};
bdgmsh.planesurface = {[1,7]; [2,8]; [3,9]; [4,10]; [5,11]; [6,12]};
for i=6+1:6+4*6
    bdgmsh.planesurface{i} = i+6;
end
bdgmsh.periodicsurface = [1,2; 3,4; 5,6];
bdgmsh.periodicsurfacelist = ... 
	{[bdgmsh.lineloop{1}, bdgmsh.lineloop{7}], ...
     [bdgmsh.lineloop{2}, bdgmsh.lineloop{8}];
     [bdgmsh.lineloop{3}, bdgmsh.lineloop{9}], ...
     [bdgmsh.lineloop{4}, bdgmsh.lineloop{10}];
     [bdgmsh.lineloop{5}, bdgmsh.lineloop{11}], ...
     [bdgmsh.lineloop{6}, bdgmsh.lineloop{12}];
     };
bdgmsh.surfaceloop = {1:numel(bdgmsh.planesurface)};
bdgmsh.volume = {1};
bdgmsh.box=[-0.5, 0.5, -0.5, 0.5, -0.5, 0.5];
end

