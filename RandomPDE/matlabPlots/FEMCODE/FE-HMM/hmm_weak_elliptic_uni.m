function [macSol, macFemspace, macMesh, macVp] = ...
    hmm_weak_elliptic_uni(macMesh, macVp, micMesh, micVp, opt) 
% HMM method for advection-diffusion multiscale problems with "weak"
% advection: -diva(a^\eps \nabla u^\eps) + b^\eps \cdot \nabla u^eps = f,
% i.e., the micro problems are purely diffusive

global global_options;
GL = global_options; % to be used in parallel workers

%% INITIALIZE
dim = size(macMesh.elem,2) - 1;

if length(macVp.elemtype) == 2
    macroSolver = @elliptic;
elseif length(macVp.elemtype) == 3 && strcmp(macVp.elemtype(3),'d')
    macroSolver = @elliptic_dg;
else
    error('chosen macro solver not implemented');
end

if length(micVp.elemtype) ~= 2
    error('chosen micro solver not implemented');
end

% only necessary to define sampling domains
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
  b0vector = zeros(NT, dim);

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% PARFOR - micro problem solutions (or updates)
  parfor k=1:NT
    
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
    
    %% CALCULATE HOMOGENIZED ADVECTION
    micVpX.b = @(x,kk) feval(micVp.b,xloc(k,:),x,kk);
    b0vector(k,:) = get_homvector_weak_adv_diff(micSolX, micFemspaceX, micMeshX, micVpX);
    
  end %% PARFOR

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  macVp.a(:,:,:,qp) = a0tensor;
  macVp.b(:,:,qp) = b0vector;
  
end
    
%% SOLVE macro
% (use the same quadrature nodes for advection --> not optimal in L2 norm)
macVp.blambda = macVp.alambda;
macVp.bweight = macVp.aweight;

[macSol, macFemspace, macMesh, macVp] = macroSolver(macMesh, macVp);
  
  
end