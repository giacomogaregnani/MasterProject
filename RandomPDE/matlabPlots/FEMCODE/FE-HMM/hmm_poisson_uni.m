function [macSol, macFemspace, macMesh, macVp, a0tensor] = ...
    hmm_poisson_uni(macMesh, macVp, micMesh, micVp, opt) 
% HMM method for purely diffusive multiscale problems

global global_options;
GL = global_options; % to be used in parallel workers

%% INITIALIZE
dim = size(macMesh.elem,2) - 1;

if length(macVp.elemtype) == 2
    macroSolver = @poisson;
elseif length(macVp.elemtype) == 3 && strcmp(macVp.elemtype(3),'d')
    macroSolver = @poisson_dg;
else
    error('chosen macro solver not implemented');
end

if length(micVp.elemtype) ~= 2
    error('chosen micro solver not implemented');
end

macFemspace = get_femspace(macMesh, macVp.elemtype);
macVp.aquad = max(2*(macFemspace.deg-1),1);
[macVp.alambda,macVp.aweight] = quadpts(dim,macVp.aquad);
macVp.fully_discrete = true;

NQ  = numel(macVp.aweight);
NT  = size(macMesh.elem,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for qp=1:NQ  % for every quadrature point
  
  xloc = get_rc(macMesh, [], 'all', [], macVp.alambda(qp,:));
  
  a0tensor = zeros(NT, dim, dim);

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% PARFOR - micro problem solutions (or updates)
  for k=1:NT
    setglobal('global_options',GL);
    fprintf('QP: %.1d / %.1d;  Elem: %.5d / %.5d\n', ...
     qp, NQ, k, NT);
    
    %% SET LOCAL VARIATIONAL PROBLEMS
    micVpX   = micVp;
    micVpX.a = @(x,kk,ll) feval(micVp.a,xloc(k,:),x,kk,ll);

    %% SOLVE
    [micSolX, micFemspaceX, micMeshX, micVpX] = poisson(micMesh, micVpX);
    
    %% CALCULATE HOMOGENIZED TENSOR
    a0tensor(k,:,:) = get_homtensor_poisson(micSolX, micFemspaceX, micMeshX, micVpX);
    
    
  end %% PARFOR

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  macVp.a(:,:,:,qp) = a0tensor;
  
end
    
%% SOLVE macro
[macSol, macFemspace, macMesh, macVp] = macroSolver(macMesh, macVp);


  
  
end