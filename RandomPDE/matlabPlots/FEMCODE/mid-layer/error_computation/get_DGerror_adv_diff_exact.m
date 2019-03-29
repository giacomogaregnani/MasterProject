function [H1norm, L2norm, H1seminorm, DGnorm] = ...
	get_DGerror_adv_diff_exact(mesh, femspace, vp, eqn, ff, f)
% calculates error in DG diffusion-advection norm compared to exact solution

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
dgsemi_diff_int = zeros(nISS,1);

for m=1:NQ_sub
    lambdam = lambda_sub(m,:); % defines integration point x_{lambdam}
    lambdam1 = lambdam(lam1(intsubsim,:));
    lambdam2 = lambdam(lam2(intsubsim,:));

    % penalization term
    dgsemi_diff_int = dgsemi_diff_int + weight_sub(m)*...
            ((evalf(mesh, subsim2elem(intsubsim,1), femspace, lambdam1, ff, 0) - ...
             evalf(mesh, subsim2elem(intsubsim,1), [], lambdam1, f, 0)) - ...
             (evalf(mesh, subsim2elem(intsubsim,2), femspace, lambdam2, ff, 0) - ...
             evalf(mesh, subsim2elem(intsubsim,2), [], lambdam2, f, 0))).^2;

end


penalization = vp.alpha./subsimdiam(intsubsim);

dgsemi_adv_int = dgsemi_diff_int.* subsimarea(intsubsim) .* eqn.flowedge(intsubsim);
dgsemi_diff_int = dgsemi_diff_int.* subsimarea(intsubsim) .* penalization;
    
nBSS = sum(isbdsubsim);
dgsemi_diff_bnd = zeros(nBSS,1);

for m=1:NQ_sub
    lambdam = lambda_sub(m,:); % defines integration point x_{lambdam}
    lambdam1 = lambdam(lam1(isbdsubsim,:));

    % penalization term
    dgsemi_diff_bnd = dgsemi_diff_bnd + weight_sub(m)*...
            ((evalf(mesh, subsim2elem(isbdsubsim,1), femspace, lambdam1, ff, 0) - ...
             evalf(mesh, subsim2elem(isbdsubsim,1), [], lambdam1, f, 0))).^2;

end


penalization = vp.alpha./subsimdiam(isbdsubsim);

dgsemi_adv_bnd = dgsemi_diff_bnd.* subsimarea(isbdsubsim) .* eqn.flowedge(isbdsubsim);
dgsemi_diff_bnd = dgsemi_diff_bnd.* subsimarea(isbdsubsim) .* penalization;

DGseminorm_diff = sqrt(sum(dgsemi_diff_int) + sum(dgsemi_diff_bnd));
DGseminorm_adv = sqrt(sum(dgsemi_adv_int) + sum(dgsemi_adv_bnd));

[H1norm, L2norm, H1seminorm] = get_H1error_exact(mesh, femspace, ff, f);

DGnorm = sqrt(eqn.scaldiff*H1seminorm^2 + eqn.scaldiff*DGseminorm_diff^2 ...
              + eqn.scaladv*L2norm^2 + DGseminorm_adv^2);


end

