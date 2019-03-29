function [bdmesh, bdgmsh]= geom3D_connected3
%GEOM_3D_CONNECTED1 Summary of this function goes here
%   Detailed explanation goes here

bdmesh = [];

%% PIPES SETTING
s = [0.1,0.1,0.1];

pipes=cell(3,1);
pipes{1} = cell(1,1);
pipes{1}{1}.coor = [-0.35,-0.35; -0.05,-0.35; -0.05,-0.05; -0.35,-0.05]; % hole
pipes{1}{1}.shift = [0,0]; % shift

pipes{2} = cell(1,1);
pipes{2}{1}.coor = [-0.35,-0.35; -0.05,-0.35; -0.05,-0.05; -0.35,-0.05]; % hole
pipes{2}{1}.shift = [0,0]; % shift

pipes{3} = cell(1,1);
pipes{3}{1}.coor = [-0.3,-0.3; -0.1,-0.3; -0.1,-0.1; -0.3,-0.1]; % hole
pipes{3}{1}.shift = [0.1,0.1]; % shift

%% GEOMETRY DESCRIPTION
verticesCube =  [-0.5,-0.5,-0.5;
 s(1),-0.5,-0.5;
-0.5, s(2),-0.5;
 s(1), s(2),-0.5;
-0.5,-0.5, s(3);
 s(1),-0.5, s(3);
-0.5, s(2), s(3);
 s(1), s(2), s(3)];

bdgmsh.node = verticesCube;
bdgmsh.line = [1,2; 3,4; 5,6; 7,8; 
	           1,3; 2,4; 5,7; 6,8;
			   1,5; 2,6; 3,7; 4,8];
bdgmsh.lineloop = {[5, 11, -7, -9]; % ACGE
	               [6, 12, -8, -10]; % BDHF	 
	               [1, 10, -3, -9]; % ABFE
	               [2, 12, -4, -11]; % CDHG
                   [1, 6, -2, -5]; %ABDC
	               [3, 8, -4, -7];} %EFGH		   
bdgmsh.planesurface = {1; 2; 3; 4; 5; 6};		
bdgmsh.periodicsurface = zeros(0,2);		   
bdgmsh.periodicsurfacelist = cell(0,2);		   

sub = [2,3; 1,3; 1,2];

for i=1:3
	for j=1:numel(pipes{i})
		% process pipes{i}{j}
		
		%% NODE
		nodeI = size(bdgmsh.node,1);
		gon = size(pipes{i}{j}.coor,1);
		
		nodes1 = repmat(s(i),[gon,3]);
		nodes1(:,sub(i,:)) = pipes{i}{j}.coor;
		
		nodes2 = repmat(0.5,[gon,3]);
		nodes2(:,sub(i,:)) = bsxfun(@plus, pipes{i}{j}.coor,pipes{i}{j}.shift);

		nodes3 = repmat(-0.5,[gon,3]);
		nodes3(:,sub(i,:)) = bsxfun(@plus, pipes{i}{j}.coor,pipes{i}{j}.shift);

		bdgmsh.node = [bdgmsh.node; nodes1; nodes2; nodes3];
		
		%% LINE
		lineI = size(bdgmsh.line,1);
		lines1 = [(nodeI+1 : nodeI+gon); (nodeI+2: nodeI+gon), nodeI+1]';
		lines2 = lines1 + gon;
		lines3 = lines1 + 2*gon;
		lines4 = [(nodeI+1: nodeI+gon); (nodeI+gon+1: nodeI+2*gon)]';
		bdgmsh.line = [bdgmsh.line; lines1; lines2; lines3; lines4];
		
		%% LINELOOP
		lineloopI = numel(bdgmsh.lineloop);
		bdgmsh.lineloop{lineloopI+1} = lineI+1 : lineI+gon;
		bdgmsh.lineloop{lineloopI+2} = lineI+gon+1 : lineI+2*gon;
		bdgmsh.lineloop{lineloopI+3} = lineI+2*gon+1 : lineI+3*gon;
		for k=1:gon
			bdgmsh.lineloop{lineloopI+3+k} = ...
				[lineI+k, lineI+3*gon+k+1, -(lineI+gon+k), -(lineI+3*gon+k)];
			if k==gon
				bdgmsh.lineloop{lineloopI+3+k}(2) = lineI + 3*gon + 1;
			end
		end
		
		%% PLANESURFACE
		holesI = numel(bdgmsh.planesurface{2*i});
		bdgmsh.planesurface{2*i}(holesI+1) = lineloopI+1; 
		bdgmsh.planesurface{2*i-1}(holesI+1) = lineloopI+3;
		
		planesurfaceI = numel(bdgmsh.planesurface);
		for k=1:gon+2
			bdgmsh.planesurface{planesurfaceI+k} = lineloopI + k + 1;
		end
		
		%% PERIODIC SURFACE
		periodicI = size(bdgmsh.periodicsurface,1);
		bdgmsh.periodicsurface(periodicI+1,1:2) = planesurfaceI + [1,2];
		bdgmsh.periodicsurfacelist{periodicI+1,1} = ...
			bdgmsh.lineloop{bdgmsh.planesurface{planesurfaceI+1}};
		bdgmsh.periodicsurfacelist{periodicI+1,2} = ...
			bdgmsh.lineloop{bdgmsh.planesurface{planesurfaceI+2}};
	end
end
		   
bdgmsh.surfaceloop = {1:numel(bdgmsh.planesurface)};
bdgmsh.volume = {1};
bdgmsh.box=[-0.5, 0.5, -0.5, 0.5, -0.5, 0.5];
end

