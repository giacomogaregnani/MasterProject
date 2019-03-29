function [stima scaldiff] = assemble_dg_diffusion_pen( mesh, femspace, vp, options )
% assembles diffusive penalization terms for DG

if nargin<4, options = struct; end

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

g = zeros(NT,dim*dim,1,NQ);
for m=1:NQ  % fourth dimension
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
            
            g(:,sub2ind([dim dim],j,k),1,m) = vpakij;

        end

    end
    
end

if isfield(options,'dgnorm') && options.dgnorm
    scaldiff = max(max(sqrt(sum(g.^2,2)),[],4),[],1);
else
    scaldiff = -1;
end

%% RESTORE VALUES FROM INTEGRATION POINTS
[rfemspace, restored] = restore(mesh, g, vp.elemtype);
clear g;

[subsim, subsim2elem, isbdsubsim] = auxstructure(mesh,'subsim','subsim2elem','isbdsubsim');
NE  = size(subsim,1);
[lambda_sub, weight_sub] = quadpts(dim-1, 2*femspace.deg);
lambda_sub(:,dim+1) = 0;
NQ_sub = numel(weight_sub);

intsubsim = ~isbdsubsim;
subsimarea = subsimplex_volume(mesh,[],[],[],subsim);
subsimdiam = subsimplex_diameter(mesh,[],[],[],subsim);

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

stima = sparse(ND,ND);

%% PENALIZATION TERM GENERATION (INTERIOR EDGES)
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
        
        % penalization term
        Aji11 = Aji11 + weight_sub(m)*...
                (evalb(mesh, subsim2elem(intsubsim,1), lambdam1, j, 0, ...
                femspace.elemtype) - ...
                0).*...
                (evalb(mesh, subsim2elem(intsubsim,1), lambdam1, i, 0, ...
                femspace.elemtype) - ...
                0);

        Aji12 = Aji12 + weight_sub(m)*...
                (evalb(mesh, subsim2elem(intsubsim,1), lambdam1, j, 0, ...
                femspace.elemtype) - ...
                0).*...
                (0 - ...
                evalb(mesh, subsim2elem(intsubsim,2), lambdam2, i, 0, ...
                femspace.elemtype));

        Aji21 = Aji21 + weight_sub(m)*...
                (0 - ...
                evalb(mesh, subsim2elem(intsubsim,2), lambdam2, j, 0, ...
                femspace.elemtype)).*...
                (evalb(mesh, subsim2elem(intsubsim,1), lambdam1, i, 0, ...
                femspace.elemtype) - ...
                0);

        Aji22 = Aji22 + weight_sub(m)*...
                (0 - ...
                evalb(mesh, subsim2elem(intsubsim,2), lambdam2, j, 0, ...
                femspace.elemtype)).*...
                (0 - ...
                evalb(mesh, subsim2elem(intsubsim,2), lambdam2, i, 0, ...
                femspace.elemtype));


    end
    
%     PROPER SACLING WITH FROBENIUS NORM OF TENSOR
    tensor1 = zeros(nISS,dim*dim);
    tensor2 = zeros(nISS,dim*dim);
    for it=1:dim^2
        tensor1(:,it) = evalf(mesh, subsim2elem(intsubsim,1), rfemspace, ...
                        lambdam1, restored(:,it,1,:));
        tensor2(:,it) = evalf(mesh, subsim2elem(intsubsim,2), rfemspace, ...
                        lambdam2, restored(:,it,1,:));
    end
    tensor_frob_norm = sqrt(sum((0.5*(tensor1 + tensor2)).^2,2));
 
    penalization = tensor_frob_norm*vp.alpha./subsimdiam(intsubsim);
    
    Aji11 = Aji11 .* subsimarea(intsubsim) .* penalization;
    Aji12 = Aji12 .* subsimarea(intsubsim) .* penalization;
    Aji21 = Aji21 .* subsimarea(intsubsim) .* penalization;
    Aji22 = Aji22 .* subsimarea(intsubsim) .* penalization;
    
        
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


%% PENALIZATION TERM GENERATION (BOUNDARY EDGES)


nBSS = sum(isbdsubsim);
for i=1:LD
for j=1:LD

    Aji11 = zeros(nBSS,1);
    
    for m=1:NQ_sub
        lambdam = lambda_sub(m,:); % defines integration point x_{lambdam}
        lambdam1 = lambdam(lam1(isbdsubsim,:));
        
        % penalization term
                
        Aji11 =  Aji11 + weight_sub(m)*...
                (evalb(mesh, subsim2elem(isbdsubsim,1), lambdam1, j, 0, ...
                femspace.elemtype) - ...
                0).*...
                (evalb(mesh, subsim2elem(isbdsubsim,1), lambdam1, i, 0, ...
                femspace.elemtype) - ...
                0);
        

    end
    
%     PROPER SACLING WITH FROBENIUS NORM OF TENSOR
    tensor1 = zeros(nBSS,dim*dim);
    for it=1:dim^2
        tensor1(:,it) = evalf(mesh, subsim2elem(isbdsubsim,1), rfemspace, ...
                        lambdam1, restored(:,it,1,:));
    end
    tensor_frob_norm = sqrt(sum(tensor1.^2,2));
    
    penalization = tensor_frob_norm*vp.alpha./subsimdiam(isbdsubsim);
    
    Aji11 = Aji11 .* subsimarea(isbdsubsim) .* penalization;
    
        
    stima = stima + sparse(get_dof(mesh, subsim2elem(isbdsubsim,1), femspace, i), ...
			get_dof(mesh, subsim2elem(isbdsubsim,1), femspace, j), Aji11, ND, ND);
    
end
end

end

