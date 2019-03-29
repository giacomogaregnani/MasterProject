function [bdmesh, bdgmsh]= geom3D_connected4
%GEOM_3D_CONNECTED1 Summary of this function goes here
%   Detailed explanation goes here

bdmesh = [];

%% PIPES SETTING
	function [s, pipes] = pore_description(x)
		s = [0.1,0.1,0.1];
		
		pipes=cell(3,1);
		pipes{1} = cell(1,1);
		pipes{1}{1}.coor = [-0.3,-0.3; -0.1,-0.3; -0.1,-0.1; -0.3,-0.1]; % hole
		pipes{1}{1}.shift = 0.1*[sin(2*pi*x(2)),cos(2*pi*x(3))]; % shift
		
		pipes{2} = cell(1,1);
		pipes{2}{1}.coor = [-0.3,-0.3; -0.1,-0.3; -0.1,-0.1; -0.3,-0.1]; % hole
		pipes{2}{1}.shift = 0.1*[sin(2*pi*x(1)),cos(2*pi*x(3))]; % shift
		
		pipes{3} = cell(1,1);
		pipes{3}{1}.coor = [-0.3,-0.3; -0.1,-0.3; -0.1,-0.1; -0.3,-0.1]; % hole
		pipes{3}{1}.shift = 0.1*[sin(2*pi*x(1)),cos(2*pi*x(2))]; % shift
	end
%[s,pipes] = pore_description([0,0,0]);

epsilon = 1/8;
llc = [0,0,0];

map = true(4,4,4);
sz = size(map);
pores = cell(sz);

bdgmsh.node = zeros(0,3);
bdgmsh.line = zeros(0,2);
bdgmsh.lineloop = cell(0,1);
bdgmsh.planesurface = cell(0,1);
bdgmsh.periodicsurface = zeros(0,2);		   
bdgmsh.periodicsurfacelist = cell(0,2);		   

for i1 = 1:sz(1)
	for i2 = 1:sz(2)
		for i3 = 1:sz(3)
			if map(i1,i2,i3)
				
				x = llc + ([i1, i2, i3]-1/2)*epsilon;
				
				[s,pipes] = pore_description(x);
				
				verticesCube =  [-0.5,-0.5,-0.5;
					s(1),-0.5,-0.5;
					-0.5, s(2),-0.5;
					s(1), s(2),-0.5;
					-0.5,-0.5, s(3);
					s(1),-0.5, s(3);
					-0.5, s(2), s(3);
					s(1), s(2), s(3)];
				verticesCube = bsxfun(@plus,verticesCube,[i1,i2,i3]-1);
				pores{i1,i2,i3}.nodeI = size(bdgmsh.node,1);
				bdgmsh.node = [bdgmsh.node; verticesCube];
				
				linesCube = pores{i1,i2,i3}.nodeI + ...
					[1,2; 3,4; 5,6; 7,8;
					1,3; 2,4; 5,7; 6,8;
					1,5; 2,6; 3,7; 4,8];
				pores{i1,i2,i3}.lineI = size(bdgmsh.line,1);
				bdgmsh.line = [bdgmsh.line; linesCube];
				
				lineloopsCube = pores{i1,i2,i3}.lineI + ...
					[5, 11, 7, 9; % ACGE
					6, 12, 8, 10; % BDHF
					1, 10, 3, 9; % ABFE
					2, 12, 4, 11; % CDHG
					1, 6, 2, 5; %ABDC
					3, 8, 4, 7]; %EFGH
				lineloopsCube(:,3:4) = -lineloopsCube(:,3:4);
				pores{i1,i2,i3}.lineloopI = numel(bdgmsh.lineloop);
				for i=1:6
					bdgmsh.lineloop{pores{i1,i2,i3}.lineloopI+i} = ...
						lineloopsCube(i,:);
				end
				
				pores{i1,i2,i3}.planesurfaceI = numel(bdgmsh.planesurface);
				for i=1:6
					bdgmsh.planesurface{pores{i1,i2,i3}.planesurfaceI+i} = ...
						pores{i1,i2,i3}.lineloopI + i;
				end				
			end
		end
	end
end

sub = [2,3; 1,3; 1,2];

for i1 = 1:sz(1)
	for i2 = 1:sz(2)
		for i3 = 1:sz(3)
			if map(i1,i2,i3)
				x = llc + ([i1, i2, i3]-1/2)*epsilon;
				[s,pipes] = pore_description(x);
				for i=1:3
					if (i==1) && (i1==sz(1)), continue; end
					if (i==2) && (i2==sz(2)), continue; end
					if (i==3) && (i3==sz(3)), continue; end
					for j=1:numel(pipes{i})
						%% NODE
						nodeI = size(bdgmsh.node,1);
						gon = size(pipes{i}{j}.coor,1);
						
						nodes1 = repmat(s(i),[gon,3]);
						nodes1(:,sub(i,:)) = pipes{i}{j}.coor;
						nodes1 = bsxfun(@plus,nodes1,[i1,i2,i3]-1);
						
						nodes2 = repmat(0.5,[gon,3]);
						nodes2(:,sub(i,:)) = bsxfun(@plus,pipes{i}{j}.coor,pipes{i}{j}.shift);
						nodes2 = bsxfun(@plus,nodes2,[i1,i2,i3]-1);
						
						bdgmsh.node = [bdgmsh.node; nodes1; nodes2]; % nodes3
						
						%% LINE
						lineI = size(bdgmsh.line,1);
						lines1 = [(nodeI+1 : nodeI+gon); (nodeI+2: nodeI+gon), nodeI+1]';
						lines2 = lines1 + gon;
						lines3 = [(nodeI+1: nodeI+gon); (nodeI+gon+1: nodeI+2*gon)]';
						bdgmsh.line = [bdgmsh.line; lines1; lines2; lines3];
						
						%% LINELOOP
						lineloopI = numel(bdgmsh.lineloop);
						bdgmsh.lineloop{lineloopI+1} = lineI+1 : lineI+gon;
						bdgmsh.lineloop{lineloopI+2} = lineI+gon+1 : lineI+2*gon;
						for k=1:gon
							bdgmsh.lineloop{lineloopI+2+k} = ...
								[lineI+k, lineI+2*gon+k+1, -(lineI+gon+k), -(lineI+2*gon+k)];
							if k==gon
								bdgmsh.lineloop{lineloopI+2+k}(2) = lineI + 2*gon + 1;
							end
						end
						
						%% PLANESURFACE						
						holesI = numel(bdgmsh.planesurface{pores{i1,i2,i3}.planesurfaceI+2*i});
						bdgmsh.planesurface{pores{i1,i2,i3}.planesurfaceI+2*i}(holesI+1) = lineloopI+1;
						j1 = i1; j2 = i2; j3 = i3;
						if i==1, j1 = j1+1; end
						if i==2, j2 = j2+1; end
						if i==3, j3 = j3+1; end
						bdgmsh.planesurface{pores{j1,j2,j3}.planesurfaceI+2*i-1}(holesI+1) = lineloopI+2;
						
						planesurfaceI = numel(bdgmsh.planesurface);
						for k=1:gon
							bdgmsh.planesurface{planesurfaceI+k} = lineloopI + k + 2;
						end
					end
				end
			end
		end
	end
end

bdgmsh.surfaceloop = {1:numel(bdgmsh.planesurface)};
bdgmsh.volume = {1};
bdgmsh.box=[-0.5, 0.5, -0.5, 0.5, -0.5, 0.5];
end

