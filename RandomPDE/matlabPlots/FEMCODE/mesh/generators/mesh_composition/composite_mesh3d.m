function mesh = composite_mesh3d(gfun, shape, M, options)
%FINE_CMESH_GENERATOR Summary of this function goes here
%   Detailed explanation goes here
if nargin<3, M=1; end
if nargin<2, shape = 'C'; end
if nargin<4, options = struct; end
if ~isfield(options,'periodic'), options.periodic = false; end
if ~islogical(shape)
	if strcmp(shape,'C')
		error('needs some updates');
		%shape = logical([1,1,0; 1,0,0; 1,1,0])';
	else
		error('Unknown shape');
	end
end

%% PREPARATION
perx=size(shape,1); peri=M*perx;
pery=size(shape,2); perj=M*pery;
perz=size(shape,3); perk=M*perz;
tol =10^-5;

subf = [-1,1,0]; subk = [2,1,3];

gmsh_options.lc = 0.15;
gmsh_options.periodic = false;
options.optimize = 3;
options.optimizenetgen = 3;

shape = repval(shape,[M,M,M]);
done =zeros(size(shape));

locmesh = cell(size(shape));

I=0;
mesh.node = [];
mesh.elem = [];
num = 0;
Ngrids = sum(sum(sum(shape)));
for i=1:peri
	for j=1:perj
		for k=1:perk
			if ~shape(i,j,k), continue; end
			
			%% DISPLAY MESSAGE
			num = num + 1;
			fprintf('Processing grid tile %d from %d\n', num, Ngrids);
			
			%% LOCALIZATION AND LOADING THE BOUNDARY
			cent = [i-1/2, j-1/2, k-1/2] / M;
			epsilon = 1/M;
			[~, bdgmsh] = gfun(cent, epsilon);
			
			lmesh = gmsh(bdgmsh, gmsh_options);
			N = size(lmesh.node,1);
			
            %% FIND BORDER
            lmesh.border = cell(3,3);
            lmesh.border{1,1} = find(sum(abs(bsxfun(@minus,lmesh.node,[-1/2,-1/2])),2) < tol);
            lmesh.border{1,2} = find(sum(abs(bsxfun(@minus,lmesh.node,[-1/2, 1/2])),2) < tol);
            lmesh.border{2,2} = find(sum(abs(bsxfun(@minus,lmesh.node,[ 1/2, 1/2])),2) < tol);
            lmesh.border{2,1} = find(sum(abs(bsxfun(@minus,lmesh.node,[ 1/2,-1/2])),2) < tol);
            for ii = 1 : 2
                lmesh.border{ii,3} = find(abs(lmesh.node(:,1) - (-1)^ii*1/2) < tol);
                [~, ind] = sort(lmesh.node(lmesh.border{ii,3},2));
                lmesh.border{ii,3} = lmesh.border{ii,3}(ind);
            end
			for ii = 1 : 2
				lmesh.border{3,ii} = find(abs(lmesh.node(:,2) - (-1)^ii*1/2) < tol);
				[~, ind] = sort(lmesh.node(lmesh.border{3,ii},1));
				lmesh.border{3,ii} = lmesh.border{3,ii}(ind);
			end
			
			%% CONTINUE AS BEFORE
			locmesh{i,j} = lmesh;
			locmesh{i,j}.node = bsxfun(@plus, locmesh{i,j}.node / M, cent);
			changed = zeros(0,2);
			for m1=1:3
				for m2=1:3
					if (m1==3) && (m2==3), continue; end
					
					i2 = i + subf(m1); j2 = j + subf(m2);
					if i2<1, i2 = peri; end
					if i2 > peri, i2 = 1; end
					if j2<1, j2 = perj; end
					if j2 > perj, j2 = 1; end
					if ~options.periodic && ((i + subf(m1) ~= i2) || ...
							(j + subf(m2) ~= j2))
						continue;
					end
					
					if ~(done(i2,j2)), continue; end
					
					
					[~,ia] = setdiff(locmesh{i,j}.border{m1,m2}, changed(:,1));
					changed = [changed; locmesh{i,j}.border{m1,m2}(ia), ...
						locmesh{i2,j2}.border{subk(m1),subk(m2)}(ia)];
				end
			end
			%% UPDATE
			unchanged = setdiff(1:N, changed(:,1))';
			done(i,j) = 1;
			
			numnew = numel(unchanged);
			newind = (I+1 : I+numnew)';
			changed = [changed; unchanged, newind];
			
			locmesh{i,j}.elem = substitute_array(...
				locmesh{i,j}.elem, ...
				changed(:,1), changed(:,2));
			for m1=1:3
				for m2=1:3
					if (m1==3) && (m2==3), continue; end
					locmesh{i,j}.border{m1,m2} = substitute_array( ...
						locmesh{i,j}.border{m1,m2}, ...
						changed(:,1), changed(:,2));
				end
			end
			
			mesh.node(newind,:) = locmesh{i,j}.node(unchanged,:);
			mesh.elem = [mesh.elem; locmesh{i,j}.elem];
			
			I = I + numnew;
		end
	end
end

%% FINAL TOUCH
mesh.box = [0,perx,0,pery];
mesh.periodic = options.periodic;

mesh = label(mesh);
mesh = fixorder(mesh);
end