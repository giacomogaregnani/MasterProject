function [bdmesh, bdgmsh] = geom3D_gen_tet_hole( tetvert )
%GEOM3D_GEN_TET_HOLE Summary of this function goes here
%   Detailed explanation goes here

A =[-0.5,-0.5,-0.5];
B =[ 0.5,-0.5,-0.5];
C =[-0.5, 0.5,-0.5];
D =[ 0.5, 0.5,-0.5];
E =[-0.5,-0.5, 0.5];
F =[ 0.5,-0.5, 0.5];
G =[-0.5, 0.5, 0.5];
H =[ 0.5, 0.5, 0.5];

bdmesh.node = [A; B; C; D; E; F; G; H; tetvert];
bdgmsh.node = bdmesh.node;
bdgmsh.line = [1,2; 3,4; 5,6; 7,8; 
	           1,3; 2,4; 5,7; 6,8;
			   1,5; 2,6; 3,7; 4,8;
			   9,10; 9,11; 9,12; 10,11; 10,12; 11,12];
bdgmsh.lineloop = {[1, 6, -2, -5]; %ABDC
	               [3, 8, -4, -7]; %EFGH
	               [1, 10, -3, -9]; % ABFE
	               [2, 12, -4, -11]; % CDHG
	               [5, 11, -7, -9]; % ACGE
	               [6, 12, -8, -10]; % BDHF
				   [13,16,-14]; 
				   [13,17,-15]; 
				   [14,18,-15]; 
				   [16,18,-17]}; %tet
bdgmsh.planesurface = {1; 2; 3; 4; 5; 6; 7; 8; 9; 10};
bdgmsh.periodicsurface = [1,2; 3,4; 5,6];
bdgmsh.periodicsurfacelist = ... 
	{[1, 6, -2, -5], [3, 8, -4, -7];
	[1, 10, -3, -9], [2, 12, -4, -11];
	[5, 11, -7, -9], [6, 12, -8, -10]};
bdgmsh.surfaceloop = {[1,2,3,4,5,6]; [7,8,9,10]};
bdgmsh.volume = {[1,2]};
bdgmsh.box=[-0.5, 0.5, -0.5, 0.5, -0.5, 0.5];

bdmesh.elem = [
	 4, 1, 2; 1, 4, 3; % ABCD: DAB + ADC
	 5, 8, 6; 8, 5, 7; % EFGH: EHF + HEG
	 1, 6, 2; 6, 1, 5; % ABFE: AFB + FAE
	 8, 3, 4; 3, 8, 7; % CDHG: HCD + CHG
	 7, 1, 3; 1, 7, 5; % ACGE: GAC + AGE
	 2, 8, 4; 8, 2, 6; % BDHF: BHD + HBF
	 9,10,12; 10,11,12; 11,9,12; 9,11,10;
	 ];

bdmesh.box=[-0.5, 0.5, -0.5, 0.5, -0.5, 0.5];

bdmesh.bdnode{1,1,1} = 1;
bdmesh.bdnode{2,1,1} = 2;
bdmesh.bdnode{1,2,1} = 3;
bdmesh.bdnode{2,2,1} = 4;
bdmesh.bdnode{1,1,2} = 5;
bdmesh.bdnode{2,1,2} = 6;
bdmesh.bdnode{1,2,2} = 7;
bdmesh.bdnode{2,2,2} = 8;

bdmesh.bdnode{3,1,1} = [1, 2];
bdmesh.bdnode{3,2,1} = [3, 4];
bdmesh.bdnode{3,1,2} = [5, 6];
bdmesh.bdnode{3,2,2} = [7, 8];

bdmesh.bdnode{1,3,1} = [1, 3];
bdmesh.bdnode{2,3,1} = [2, 4];
bdmesh.bdnode{1,3,2} = [5, 7];
bdmesh.bdnode{2,3,2} = [6, 8];

bdmesh.bdnode{1,1,3} = [1, 5];
bdmesh.bdnode{2,1,3} = [2, 6];
bdmesh.bdnode{1,2,3} = [3, 7];
bdmesh.bdnode{2,2,3} = [4, 8];

bdmesh.bdnode{3,3,1} = [1, 2, 3, 4];
bdmesh.bdnode{3,3,2} = [5, 6, 7, 8];

bdmesh.bdnode{3,1,3} = [1, 2, 5, 6];
bdmesh.bdnode{3,2,3} = [3, 4, 7, 8];

bdmesh.bdnode{1,3,3} = [1, 3, 5, 7];
bdmesh.bdnode{2,3,3} = [2, 4, 6, 8];

bdmesh.bdelem{3,3,1} = [1, 2];
bdmesh.bdelem{3,3,2} = [3, 4];

bdmesh.bdelem{3,1,3} = [5, 6];
bdmesh.bdelem{3,2,3} = [7, 8];

bdmesh.bdelem{1,3,3} = [9, 10];
bdmesh.bdelem{2,3,3} = [11, 12];

end

