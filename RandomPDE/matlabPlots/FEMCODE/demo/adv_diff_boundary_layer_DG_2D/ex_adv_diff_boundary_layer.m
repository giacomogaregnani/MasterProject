function xx = ex_adv_diff_boundary_layer(x,scal)
%% EXACT SOLUTION FOR A 2D ADVECTION-DIFFUSION PROBLEM WITH BOUNDARY LAYER

% ONLY u(x) POSSIBLE
u_add = @(x,y) x.*y.*(1+(exp(-1/scal)-exp(-(1-x).*(1-y)/scal))/(1-exp(-1/scal)));
xx = u_add(x(:,1),x(:,2));
    
end

