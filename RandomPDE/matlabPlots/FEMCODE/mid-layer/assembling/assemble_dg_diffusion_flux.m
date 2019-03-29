function stima = assemble_dg_diffusion_flux( mesh, femspace, vp )
% assembles diffusive flux terms for DG

%% INITIALIZATION
dim = size(mesh.elem,2) - 1;
NT = size(mesh.elem,1);
ND = femspace.ndof;
LD = femspace.ldof;
lambda = vp.alambda;
weight = vp.aweight;
NQ = numel(weight);

%% INIT
if ~isnumeric(vp.a)
	xloc = cell([1, NQ]);
	for m=1:numel(weight)
		xloc{m} = get_rc(mesh,[],'all',[],lambda(m,:));
	end
end

g = zeros(NT,dim,LD,NQ);
for m=1:NQ  % fourth dimension
    for i=1:LD  % third dimension
        for j=1:dim  % second dimension
            
            for k=1:dim  % summing
                
                if isnumeric(vp.a)
                    if vp.fully_discrete % fully discrete tensor
                        vpakij = vp.a(:,j,k,m);
                    else
                        vpakij = vp.a(j,k);
                    end
                else
                    vpakij = vp.a(get_rc(mesh,[],'all',[],lambda(m,:)),j,k);
                end
                
                % values of a*grad v^H
                g(:,j,i,m) = g(:,j,i,m) + vpakij .* ...
                        evalb(mesh, 'all', lambda(m,:), i, k, femspace.elemtype);	
                    
                
            end
            
        end
        
    end
    
end

%% RESTORE VALUES FROM INTEGRATION POINTS
[rfemspace, restored] = restore(mesh, g, vp.elemtype);
clear g;

% quadrature formula on subsims (last weight is set tozero)
[normals, subsim, subsim2elem] = get_normals(mesh);
NE  = size(subsim,1);
[lambda_sub, weight_sub] = quadpts(dim-1, rfemspace.deg+femspace.deg); % as integrate product rfemspace.deg*femspace.deg
lambda_sub(:,dim+1) = 0;
NQ_sub = numel(weight_sub);

isbdsubsim = auxstructure(mesh,'isbdsubsim');
intsubsim = ~isbdsubsim;
subsimarea = subsimplex_volume(mesh,[],[],[],subsim);

%% Arrays lam1 and lam2 - computation
% for details: see poisson_jump_residual.m
lam1 = repmat(dim+1,[NE,dim+1]);
lam2 = repmat(dim+1,[NE,dim+1]);
for i=1:dim+1
	for j=1:dim
		ix1 = subsim(:,j) == mesh.elem(subsim2elem(:,1),i);
		ix2 = subsim(:,j) == mesh.elem(subsim2elem(:,2),i);
		lam1(ix1,i) = j;
		lam2(ix2,i) = j;
	end
end
clear ix1 ix2;

stima=sparse(ND, ND);

%% FLUX TERM GENERATION (INTERIOR EDGES)
nISS = sum(intsubsim);

for i=1:LD
for j=1:LD

    Aji11 = zeros(nISS,1);
    Aji12 = zeros(nISS,1);
    Aji21 = zeros(nISS,1);
    Aji22 = zeros(nISS,1);
    
    for m=1:NQ_sub
        lambdam = lambda_sub(m,:); % defines integration point x_{lambdam}
        lambdam1 = lambdam(lam1(intsubsim,:));
        lambdam2 = lambdam(lam2(intsubsim,:));
        
        % flux term
        for k=1:dim
                
            Aji11 = Aji11 + weight_sub(m)*bsxfun(@times, ...
                    0.5*(evalf(mesh, subsim2elem(intsubsim,1), rfemspace, ...
                    lambdam1, restored(:,k,j,:)) + ...
                    0), ...
                    normals(intsubsim,k)).*...
                    (evalb(mesh, subsim2elem(intsubsim,1), lambdam1, i, 0, ...
                    femspace.elemtype) - ...
                    0);
    
            Aji12 = Aji12 + weight_sub(m)*bsxfun(@times, ...
                    0.5*(evalf(mesh, subsim2elem(intsubsim,1), rfemspace, ...
                    lambdam1, restored(:,k,j,:)) + ...
                    0), ...
                    normals(intsubsim,k)).*...
                    (0 - ...
                    evalb(mesh, subsim2elem(intsubsim,2), lambdam2, i, 0, ...
                    femspace.elemtype));
            
            Aji21 = Aji21 + weight_sub(m)*bsxfun(@times, ...
                    0.5*(0 + ...
                    evalf(mesh, subsim2elem(intsubsim,2), rfemspace, ...
                    lambdam2, restored(:,k,j,:))), ...
                    normals(intsubsim,k)).*...
                    (evalb(mesh, subsim2elem(intsubsim,1), lambdam1, i, 0, ...
                    femspace.elemtype) - ...
                    0);
                
            Aji22 = Aji22 + weight_sub(m)*bsxfun(@times, ...
                    0.5*(0 + ...
                    evalf(mesh, subsim2elem(intsubsim,2), rfemspace, ...
                    lambdam2, restored(:,k,j,:))), ...
                    normals(intsubsim,k)).*...
                    (0 - ...
                    evalb(mesh, subsim2elem(intsubsim,2), lambdam2, i, 0, ...
                    femspace.elemtype));
        end

    end
    
    Aji11 = Aji11 .* subsimarea(intsubsim);
    Aji12 = Aji12 .* subsimarea(intsubsim);
    Aji21 = Aji21 .* subsimarea(intsubsim);
    Aji22 = Aji22 .* subsimarea(intsubsim);
    
        
    stima = stima + sparse(get_dof(mesh, subsim2elem(intsubsim,1), femspace, i), ...
			get_dof(mesh, subsim2elem(intsubsim,1), femspace, j), Aji11, ND, ND);
    
    stima = stima + sparse(get_dof(mesh, subsim2elem(intsubsim,2), femspace, i), ...
			get_dof(mesh, subsim2elem(intsubsim,1), femspace, j), Aji12, ND, ND);
        
    stima = stima + sparse(get_dof(mesh, subsim2elem(intsubsim,1), femspace, i), ...
			get_dof(mesh, subsim2elem(intsubsim,2), femspace, j), Aji21, ND, ND);
        
    stima = stima + sparse(get_dof(mesh, subsim2elem(intsubsim,2), femspace, i), ...
			get_dof(mesh, subsim2elem(intsubsim,2), femspace, j), Aji22, ND, ND);
end
end

%% FLUX TERM GENERATION (BOUNDARY EDGES)

% stima_temp=sparse(ND, ND);
nBSS = sum(isbdsubsim);
for i=1:LD
for j=1:LD

    Aji11 = zeros(nBSS,1);
    
    for m=1:NQ_sub
        lambdam = lambda_sub(m,:); % defines integration point x_{lambdam}
        lambdam1 = lambdam(lam1(isbdsubsim,:));
        
        % flux term
        for k=1:dim
                
            Aji11 = Aji11 + weight_sub(m)*bsxfun(@times, ...
                    evalf(mesh, subsim2elem(isbdsubsim,1), rfemspace, ...
                    lambdam1, restored(:,k,j,:)), ...
                    normals(isbdsubsim,k)).*...
                    evalb(mesh, subsim2elem(isbdsubsim,1), lambdam1, i, 0, ...
                    femspace.elemtype);
        end

    end
    
    Aji11 = Aji11 .* subsimarea(isbdsubsim);
    
        
    stima = stima + sparse(get_dof(mesh, subsim2elem(isbdsubsim,1), femspace, i), ...
			get_dof(mesh, subsim2elem(isbdsubsim,1), femspace, j), Aji11, ND, ND);
    
end
end


end

