function [stima scaladv flowedge] = assemble_dg_advection_pen( mesh, femspace, vp, options )
% Assembles DG penalization matrix for advection
% implemented only for P1-DG

if nargin<4, options = struct; end

%% INITIALIZATION
dim = size(mesh.elem,2) - 1;
NT = size(mesh.elem,1);
ND = femspace.ndof;
LD = femspace.ldof;
lambda = vp.alambda;
weight = vp.aweight;
NQ = numel(weight);

if LD ~= dim+1 || NQ ~= 1
    error('only p1 implemented for dg advection');
end

%% INIT
if ~isnumeric(vp.a)
	xloc = cell([1, NQ]);
	for m=1:numel(weight)
		xloc{m} = get_rc(mesh,[],'all',[],lambda(m,:));
	end
end


belem = zeros(NT,dim);

% formulated only for P1-DGFEM --> only 1 quadrature point, m=1
if numel(weight)>1, error('DG for advection only with P1 elements'); end
for j=1:dim  % second dimension
    
    if isnumeric(vp.b)
        if vp.fully_discrete
            vpb = vp.b(:,j,1);
        else
            vpb = vp.b(j);
        end
    else
        vpb = vp.b(xloc{1}, j);
    end
    
    belem(:,j) = vpb;

end
        

if options.dgnorm
    scaladv = max(sqrt(sum(belem.^2,2)),[],1);
else
    scaladv = -1;
end

% quadrature formula on subsims (last weight is set tozero)
[normals, subsim, subsim2elem] = get_normals(mesh);
NE  = size(subsim,1);
[lambda_sub, weight_sub] = quadpts(dim-1, 2*femspace.deg); % as integrate product femspace.deg*femspace.deg
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

%% ADVECTION PENALIZATION TERM GENERATION (INTERIOR EDGES)
nISS = sum(intsubsim);

bav = 0.5*(belem(subsim2elem(:,1),:) + belem(subsim2elem(:,2),:)).*repmat(intsubsim,1,dim) + ...
    belem(subsim2elem(:,1),:).*repmat(isbdsubsim,1,dim);
flowedge = sum(bav.*normals,2);
inflowelem = subsim2elem(:,1);
inflowelem(flowedge > 0) = subsim2elem(flowedge > 0,2).*uint32(intsubsim(flowedge > 0));
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


        Aji11 = Aji11 + weight_sub(m)*flowedge(intsubsim).*...
            (evalb(mesh, subsim2elem(intsubsim,1), lambdam1, j, 0, femspace.elemtype) - ...
            0).*...
            evalb(mesh, subsim2elem(intsubsim,1), lambdam1, i, 0, femspace.elemtype).*...
            (subsim2elem(intsubsim,1) == inflowelem(intsubsim));


        Aji12 = Aji12 + weight_sub(m)*flowedge(intsubsim).*...
            (evalb(mesh, subsim2elem(intsubsim,1), lambdam1, j, 0, femspace.elemtype) - ...
            0).*...
            evalb(mesh, subsim2elem(intsubsim,2), lambdam2, i, 0, femspace.elemtype).*...
            (subsim2elem(intsubsim,2) == inflowelem(intsubsim));

        Aji21 = Aji21 + weight_sub(m)*flowedge(intsubsim).*...
            (0 - ...
            evalb(mesh, subsim2elem(intsubsim,2), lambdam2, j, 0, femspace.elemtype)).*...
            evalb(mesh, subsim2elem(intsubsim,1), lambdam1, i, 0, femspace.elemtype).*...
            (subsim2elem(intsubsim,1) == inflowelem(intsubsim));

        Aji22 = Aji22 + weight_sub(m)*flowedge(intsubsim).*...
            (0 - ...
            evalb(mesh, subsim2elem(intsubsim,2), lambdam2, j, 0, femspace.elemtype)).*...
            evalb(mesh, subsim2elem(intsubsim,2), lambdam2, i, 0, femspace.elemtype).*...
            (subsim2elem(intsubsim,2) == inflowelem(intsubsim));

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
        
                
        Aji11 = Aji11 + weight_sub(m)*flowedge(isbdsubsim).*...
        (evalb(mesh, subsim2elem(isbdsubsim,1), lambdam1, j, 0, femspace.elemtype) - ...
        0).*...
        evalb(mesh, subsim2elem(isbdsubsim,1), lambdam1, i, 0, femspace.elemtype).*...
        (subsim2elem(isbdsubsim,1) == inflowelem(isbdsubsim));

    end
    
    Aji11 = Aji11 .* subsimarea(isbdsubsim);
    
        
    stima = stima + sparse(get_dof(mesh, subsim2elem(isbdsubsim,1), femspace, i), ...
			get_dof(mesh, subsim2elem(isbdsubsim,1), femspace, j), Aji11, ND, ND);
    
end
end


end

