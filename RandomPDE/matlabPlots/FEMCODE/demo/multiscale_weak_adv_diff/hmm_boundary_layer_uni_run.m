%% USER INPUT

micElem = 'p1';
Nmic = 4;
macElem  = 'p1d';
Nmac = 8;
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

expName  = 'BoundaryLayerMS2D';
dim = 2;

%% micro problem definition
clear micTensor micVp micMesh;
micTensor = str2func(['aeps_' expName]);
micVector = str2func(['beps_' expName]);
if strcmp(coupling,'periodic')
    micVp = struct('elemtype', micElem, 'a', micTensor, 'b', micVector);
    micMesh = structured_mesh(ratio*repmat([0,1],[1,dim]),Nmic,struct('periodic',true, 'centre', true));
elseif strcmp(coupling,'dirichlet');
    micVp = struct('elemtype', micElem, 'bc', 0, 'a', micTensor, 'b', micVector);
    micMesh = structured_mesh(ratio*repmat([0,1],[1,dim]),Nmic,struct('centre',true));
else
    error('unknown multiscale coupling');
end
micVp.f = permute([zeros(dim,1), -eye(dim)],[3,2,1]);
micMesh.bdflag = 'dirichlet'; % applies to all the boundary = zero boundary

%% macro problem definition

clear macRhs macVp macMesh macSol macFemspace;
macRhs = str2func(['rhs_' expName]);
macVp = struct('elemtype', macElem, 'bc', 0, 'f', macRhs);

% specify DG parameters
if length(macElem) == 3 && strcmp(macElem(3),'d')
    if ~exist('diffpart','var'),   diffpart = 'iipg'; ...
        warning('IIPG for diffusion, as not specified'); end
    if ~exist('alpha','var'),   alpha = 10; ...
        warning('default penalization parameter, as not specified'); end
    macVp.diffpart = diffpart;
    macVp.alpha = alpha;
end

macMesh = structured_mesh(repmat([0,1],[1,dim]),Nmac,struct('centre',true));
macMesh.bdflag = 'dirichlet';

%% run HMM
options = struct();

[macSol, macFemspace, macMesh, macVp] = ...
          hmm_weak_elliptic_uni(macMesh, macVp, micMesh, micVp, options);

if length(macElem) == 2      
    simpplot_sol(macMesh, macSol);
elseif length(macElem) == 3 && strcmp(macElem(3),'d')
    simpplot_sol_dg(macMesh, macSol, macFemspace);
end