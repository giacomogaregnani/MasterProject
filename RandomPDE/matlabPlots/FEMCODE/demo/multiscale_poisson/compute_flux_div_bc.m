function flux = compute_flux_div_bc(idx)
micElem = 'p1';
Nmic = 64;
macElem  = 'p1';
Nmac = 64;
ratio = 1; % sampling size: delta/epsilon
coupling = 'periodic';
% coupling = 'dirichlet';

%% parameters for DG macro solver (not mandatory)
% diffpart = 'sipg'; % SYMMETRIC INTERIOR PENALTY
% diffpart = 'nipg'; % NON-SYMMETRIC INTERIOR PENALTY
% diffpart = 'iipg'; % INCOMPLETE INTERIOR PENALTY

% penalization parameter
% alpha = 10;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MODIFICATION BELOW THIS LINE ONLY TO CHANGE PROBLEM'S DEFINITION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

expName  = 'PoissonMS2D';
dim = 2;

%% micro problem definition
clear micTensor micVp micMesh;
micTensor = str2func(['aeps_' expName]);
if strcmp(coupling,'periodic')
    micVp = struct('elemtype', micElem, 'a', micTensor);
    micMesh = structured_mesh(ratio*repmat([0,1],[1,dim]),Nmic,struct('periodic',true, 'centre', false)); % if true refine at centre
elseif strcmp(coupling,'dirichlet');
    micVp = struct('elemtype', micElem, 'bc', 0, 'a', micTensor);
    micMesh = structured_mesh(ratio*repmat([0,1],[1,dim]),Nmic,struct('centre',true));
else
    error('unknown multiscale coupling');
end
micVp.f = permute([zeros(dim,1), -eye(dim)],[3,2,1]);
micMesh.bdflag = 'dirichlet'; % applies to all the boundary = zero boundary

%% macro problem definition

clear macRhs macVp macMesh macSol macFemspace;
macRhs = str2func(['rhs_' expName]);
macVp = struct('elemtype', macElem, 'bc', @(x) dirichlet(x,idx), 'f', macRhs); % bc -> boundary condition ?

% specify DG parameters
if length(macElem) == 3 && strcmp(macElem(3),'d')
    if ~exist('diffpart','var'),   diffpart = 'sipg'; ...
        warning('SIPG for diffusion, as not specified'); end
    if ~exist('alpha','var'),   alpha = 10; ...
        warning('default penalization parameter, as not specified'); end
    macVp.diffpart = diffpart;
    macVp.alpha = alpha;
end

macMesh = structured_mesh(repmat([0,1],[1,dim]),Nmac,struct('centre',false)); % if true refine at center
macMesh.bdflag = 'dirichlet';

%% run HMM
options = struct();

[macSol, macFemspace, macMesh, macVp, a0tensor] = ...
          hmm_poisson_uni(macMesh, macVp, micMesh, micVp, options);
      
h = 1/Nmac;
A = assemble_stima(macMesh, macFemspace, macVp);
A_restricted = A(macFemspace.bd_dof,:);
residual = A_restricted*macSol;
rhs = residual;
flux = rhs/h;
n_bd_nodes = length(flux);
indexing = [1:n_bd_nodes/4, n_bd_nodes/4+1:2:n_bd_nodes*3/4-1, n_bd_nodes:-1:n_bd_nodes*3/4, n_bd_nodes*3/4-2:-2:n_bd_nodes/4+2];
flux = flux(indexing);
flux = [flux; flux(1)];

% if length(macElem) == 2      
%     simpplot_sol(macMesh, macSol);
% elseif length(macElem) == 3 && strcmp(macElem(3),'d')
%     simpplot_sol_dg(macMesh, macSol, macFemspace);
% end