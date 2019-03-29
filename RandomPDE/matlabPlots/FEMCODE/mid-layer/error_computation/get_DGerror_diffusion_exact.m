function [H1norm, L2norm, H1seminorm, DGnorm, DGseminorm] = ...
	get_DGerror_diffusion_exact(mesh, femspace, vp, ff, f)
% calculates error in DG diffusion norm compared to exact solution

dim = size(mesh.elem,2) - 1;
deg = 2*femspace.deg;

[subsim, subsim2elem, isbdsubsim] = auxstructure(mesh,'subsim','subsim2elem','isbdsubsim');
NE  = size(subsim,1);
[lambda_sub, weight_sub] = quadpts(dim-1, deg);
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


%% CALCULATE NORM OF ERROR (INTERIOR EDGES)
nISS = sum(intsubsim);
dgsemi_int = zeros(nISS,1);


for m=1:NQ_sub
    lambdam = lambda_sub(m,:); % defines integration point x_{lambdam}
    lambdam1 = lambdam(lam1(intsubsim,:));
    lambdam2 = lambdam(lam2(intsubsim,:));

    % penalization term
    dgsemi_int = dgsemi_int + weight_sub(m)*...
            ((evalf(mesh, subsim2elem(intsubsim,1), femspace, lambdam1, ff, 0) - ...
             evalf(mesh, subsim2elem(intsubsim,1), [], lambdam1, f, 0)) - ...
             (evalf(mesh, subsim2elem(intsubsim,2), femspace, lambdam2, ff, 0) - ...
             evalf(mesh, subsim2elem(intsubsim,2), [], lambdam2, f, 0))).^2;

end


penalization = vp.alpha./subsimdiam(intsubsim);


dgsemi_int = dgsemi_int.* subsimarea(intsubsim) .* penalization;
    
nBSS = sum(isbdsubsim);
dgsemi_bnd = zeros(nBSS,1);

for m=1:NQ_sub
    lambdam = lambda_sub(m,:); % defines integration point x_{lambdam}
    lambdam1 = lambdam(lam1(isbdsubsim,:));

    % penalization term
    dgsemi_bnd = dgsemi_bnd + weight_sub(m)*...
            ((evalf(mesh, subsim2elem(isbdsubsim,1), femspace, lambdam1, ff, 0) - ...
             evalf(mesh, subsim2elem(isbdsubsim,1), [], lambdam1, f, 0))).^2;

end


penalization = vp.alpha./subsimdiam(isbdsubsim);


dgsemi_bnd = dgsemi_bnd.* subsimarea(isbdsubsim) .* penalization;

DGseminorm = sqrt(sum(dgsemi_int) + sum(dgsemi_bnd));

[H1norm, L2norm, H1seminorm] = get_H1error_exact(mesh, femspace, ff, f);

DGnorm = sqrt(H1seminorm^2 + DGseminorm^2);


end

