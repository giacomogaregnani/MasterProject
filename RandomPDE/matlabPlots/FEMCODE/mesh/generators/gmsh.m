function mesh = gmsh(bdgmsh, options)
%GMSH Summary of this function goes here
%   Detailed explanation goes here

global global_options;
%% PARAMETER READING
d = size(bdgmsh.node,2);
if nargin<2, options = struct; end
if ~exist('global_options','var') || isempty(global_options)
	global_options.gmsh = '~/gmsh/bin/./gmsh';
	global_options.auxdata = '~/repos/experiments/auxdata/';
	global_options.agmg = true;
end
if ~isfield(options,'lc'), options.lc =0.2; end
if ~isfield(options,'periodic'), options.periodic = true; end
if ~isfield(options,'dim'), options.dim = d; end
if ~isfield(options,'clean'), options.clean = true; end
if ~isfield(options,'optimize'), options.optimize = 1; end
if ~isfield(options,'verbose'), options.verbose = false; end
if ~isfield(options,'optimizenetgen'), options.optimizenetgen = 1; end

%% FILES
while true
    if ~isfield(options,'outputname') || isempty(options.outputname),
        options.outputname = ['mesh' char(java.util.UUID.randomUUID())];
    end
    geofile = [global_options.auxdata options.outputname '.geo'];
    mshfile = [global_options.auxdata options.outputname '.msh'];
    if ~exist(geofile,'file') && ~exist(mshfile,'file')
        id = fopen(geofile,'w');
        if id > 3
            break;
        else
            pause(0.001);
        end
    else
        options.outputname = '';
    end
end

%% CONSOLE COMMAND
consolecommand = [global_options.gmsh ' ' geofile ... % INPUT GEO FILE
    ' -' num2str(options.dim) ... % MESHING
    ' -o ' mshfile]; % OUTPUT MSH FILE

%% MESH SIZE
%fprintf('Opening: %s\n',geofile);
fprintf(id,'lc = %f;\n',options.lc);

%% POINTS
for i=1:size(bdgmsh.node,1)
	fprintf(id,'Point(%d) = {',i);
	for j=1:d
		fprintf(id,'%f,', bdgmsh.node(i,j));
	end
	if (d == 2)
		fprintf(id,'0,');
	end
		fprintf(id,'lc};\n');	
end
fprintf(id,'\n');

%% LINES
for i=1:size(bdgmsh.line,1)
	fprintf(id,'Line(%d) = {%d, %d};\n', ...
		i,bdgmsh.line(i,1), bdgmsh.line(i,2));	
end
fprintf(id,'\n');

%% PERIODIC LINES
if options.periodic && isfield(bdgmsh,'periodicline')
	for i=1:size(bdgmsh.periodicline,1);
		fprintf(id,'Periodic Line {%d} = {%d};\n', ...
			bdgmsh.periodicline(i,1), bdgmsh.periodicline(i,2));
	end
	fprintf(id,'\n');
end

%% LINE LOOPS
if isfield(bdgmsh,'lineloop')
    for i=1:numel(bdgmsh.lineloop)
    	fprintf(id,'Line Loop(%d) = {%d',i, bdgmsh.lineloop{i}(1));
        for j=2:numel(bdgmsh.lineloop{i})
    		fprintf(id,', %d', bdgmsh.lineloop{i}(j));
        end
    	fprintf(id,'};\n');	
    end
    fprintf(id,'\n');
end

%% PLANE SURFACES
if isfield(bdgmsh,'planesurface')
    for i=1:numel(bdgmsh.planesurface)
    	fprintf(id,'Plane Surface(%d) = {%d',i, bdgmsh.planesurface{i}(1));
        for j=2:numel(bdgmsh.planesurface{i})
    		fprintf(id,', %d', bdgmsh.planesurface{i}(j));
        end
    	fprintf(id,'};\n');	
    end
    fprintf(id,'\n');
end

%% PERIODIC SURFACE
if options.periodic && isfield(bdgmsh,'periodicsurface')
	for i=1:size(bdgmsh.periodicsurface,1)
		fprintf(id,'Periodic Surface %d {%d', ...
			bdgmsh.periodicsurface(i,1), ...
			bdgmsh.periodicsurfacelist{i,1}(1));
		for j=2:numel(bdgmsh.periodicsurfacelist{i,1})
			fprintf(id,', %d', ...
				bdgmsh.periodicsurfacelist{i,1}(j));
		end
		fprintf(id,'} = %d {%d', ...
			bdgmsh.periodicsurface(i,2), ...
			bdgmsh.periodicsurfacelist{i,2}(1));
		for j=2:numel(bdgmsh.periodicsurfacelist{i,2})
			fprintf(id,', %d', ...
				bdgmsh.periodicsurfacelist{i,2}(j));
		end
		fprintf(id,'};\n');
	end
	fprintf(id,'\n');	
end

%% SURFACE LOOP
if isfield(bdgmsh,'surfaceloop')
	for i=1:numel(bdgmsh.surfaceloop)
		fprintf(id,'Surface Loop(%d) = {%d',i, bdgmsh.surfaceloop{i}(1));
		for j=2:numel(bdgmsh.surfaceloop{i})
			fprintf(id,', %d', bdgmsh.surfaceloop{i}(j));
		end
		fprintf(id,'};\n');
	end
	fprintf(id,'\n');
end

%% VOLUME
if isfield(bdgmsh,'volume')
	for i=1:numel(bdgmsh.volume)
		fprintf(id,'Volume(%d) = {%d',i, bdgmsh.volume{i}(1));
		for j=2:numel(bdgmsh.volume{i})
			fprintf(id,', %d', bdgmsh.volume{i}(j));
		end
		fprintf(id,'};\n');
	end
	fprintf(id,'\n');
end

%% POST PROCESSING - OPTIMIZE MESH
if options.optimize>0
	fprintf(id,'Mesh.Optimize = %d;\n',options.optimize);
end
if options.optimizenetgen > 0
	fprintf(id,'Mesh.OptimizeNetgen = %d;\n', options.optimizenetgen);
end

fclose(id); 

%% MESH USING GMESH
[~,~] = system(consolecommand);

%% READ .MSH FILE
mesh = read_msh(mshfile, d, options);

%% FIND PERIODIC NODES BY BRUTE FORCE
if options.periodic
	mesh = periodize(mesh, bdgmsh.box);
end

mesh = label(mesh);
mesh = fixorder(mesh);

%% DELETE AUXILIARY FILES
if options.clean
	delete(geofile);
	delete(mshfile);
end
end

