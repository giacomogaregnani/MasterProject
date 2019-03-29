function mesh = read_msh(mshfile, d, options)
%READ_MSH Summary of this function goes here
%   Detailed explanation goes here

mesh = struct;
id = fopen(mshfile, 'r');

fgetl(id); %MeshFormat
fscanf(id, '%f %f %f\n',[1,3]); % FORMAT
a=fgetl(id); %EndMeshFormat
while ~strcmp(a,'$Nodes')
    a=fgetl(id); %Nodes
end
NP = fscanf(id, '%d', [1 1]);

%% SCAN NODES
mesh.node = fscanf(id, '%f', [4 NP]);
mesh.node = mesh.node';
if (d==2)
	mesh.node(:,[1,4]) = [];
elseif (d==3)
	mesh.node(:,1) = [];
end

fgetl(id); % empty
fgetl(id); % EndNodes
fgetl(id); % Elements
NE = fscanf(id, '%d', [1 1]);

%% SCAN ELEMENTS
NT = 0;
mesh.elem = zeros(NE,options.dim + 1);
for i=1:NE
	if options.verbose && mod(i,1000) == 0;
		fprintf('Procesing: %d / %d\n',i,NE);
	end
	var = fscanf(id, '%d', [3,1]);
	type = var(2);
	numtag = var(3);
	fscanf(id, '%d', [numtag,1]);
	switch type 
		case 15
			fscanf(id, '%d', [1,1]);
		case 1
            if (options.dim == 1)
				NT = NT + 1;
				mesh.elem(NT,1:3) = fscanf(id, '%d', [1,2]);
            else
                fscanf(id, '%d', [1,2]);
            end
		case 2
			if (options.dim == 2)
				NT = NT + 1;
				mesh.elem(NT,1:3) = fscanf(id, '%d', [1,3]);
			else
				fscanf(id, '%d', [1,3]);
			end
		case 4
			if (options.dim == 3)
				NT = NT + 1;
				mesh.elem(NT,1:4) = fscanf(id, '%d', [1,4]);
			else
				fscanf(id, '%d', [1,4]);
			end
	end
end
fclose(id);
mesh.elem = mesh.elem(1:NT,:);

end

