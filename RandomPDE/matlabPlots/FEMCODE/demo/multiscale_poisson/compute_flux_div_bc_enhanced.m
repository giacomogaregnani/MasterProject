function flux = compute_flux_div_bc_enhanced(idx)
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
% flux on edge y = 0
number_internal_nodes = Nmac - 1;
[index_y0,~]=find(macMesh.node(macFemspace.bd_dof,2)==0);
flux_at_first_node = ([reshape(macVp.a(2,1,:),1,2); reshape(macVp.a(2,2,:),1,2)]*grad_dirichlet(macMesh.node(macFemspace.bd_dof(index_y0(1)),:),idx))'*[0; -1];
flux_at_last_node = ([reshape(macVp.a(Nmac*2-1,1,:),1,2); reshape(macVp.a(Nmac*2-1,2,:),1,2)]*grad_dirichlet(macMesh.node(macFemspace.bd_dof(index_y0(end)),:),idx))'*[0; -1];
A = assemble_stima(macMesh, macFemspace, macVp);
A_restricted = A(macFemspace.bd_dof(index_y0(2:end-1)),:);
residual = A_restricted*macSol;
rhs = residual - (h/6)*[flux_at_first_node; zeros(number_internal_nodes-2,1); flux_at_last_node];
T = (h/6)*(diag(4*ones(number_internal_nodes,1)) + diag(ones(number_internal_nodes-1,1),1) + diag(ones(number_internal_nodes-1,1),-1));
flux = T\rhs;
f1 = [flux_at_first_node; flux; flux_at_last_node];

% flux on edge x = 1
[index_x1,~]=find(macMesh.node(macFemspace.bd_dof,1)==1);
flux_at_first_node = ([reshape(macVp.a(Nmac*2-1,1,:),1,2); reshape(macVp.a(Nmac*2-1,2,:),1,2)]*grad_dirichlet(macMesh.node(macFemspace.bd_dof(index_x1(1)),:),idx))'*[1; 0];
flux_at_last_node = ([reshape(macVp.a(2*Nmac*Nmac,1,:),1,2); reshape(macVp.a(2*Nmac*Nmac,2,:),1,2)]*grad_dirichlet(macMesh.node(macFemspace.bd_dof(index_x1(end)),:),idx))'*[1; 0];
A_restricted = A(macFemspace.bd_dof(index_x1(2:end-1)),:);
residual = A_restricted*macSol;
rhs =  residual - (h/6)*[flux_at_first_node; zeros(number_internal_nodes-2,1); flux_at_last_node];
flux = T\rhs;
f2 = [flux_at_first_node; flux; flux_at_last_node];

% flux on edge y = 1
[index_y1,~]=find(macMesh.node(macFemspace.bd_dof,2)==1);
flux_at_first_node = ([reshape(macVp.a(2*Nmac*(Nmac-1)+2,1,:),1,2); reshape(macVp.a(2*Nmac*(Nmac-1)+2,2,:),1,2)]*grad_dirichlet(macMesh.node(macFemspace.bd_dof(index_y1(1)),:),idx))'*[0; 1];
flux_at_last_node = ([reshape(macVp.a(2*Nmac*Nmac-1,1,:),1,2); reshape(macVp.a(2*Nmac*Nmac-1,2,:),1,2)]*grad_dirichlet(macMesh.node(macFemspace.bd_dof(index_y1(end)),:),idx))'*[0; 1];
A_restricted = A(macFemspace.bd_dof(index_y1(2:end-1)),:);
residual = A_restricted*macSol;
rhs = residual - (h/6)*[flux_at_first_node; zeros(number_internal_nodes-2,1); flux_at_last_node];
flux = T\rhs;
f3 = [flux_at_first_node; flux; flux_at_last_node];

% flux on edge x = 0
[index_x0,~]=find(macMesh.node(macFemspace.bd_dof,1)==0);
flux_at_first_node = ([reshape(macVp.a(1,1,:),1,2); reshape(macVp.a(1,2,:),1,2)]*grad_dirichlet(macMesh.node(macFemspace.bd_dof(index_x0(1)),:),idx))'*[-1; 0];
flux_at_last_node = ([reshape(macVp.a(2*Nmac*(Nmac-1)+2,1,:),1,2); reshape(macVp.a(2*Nmac*(Nmac-1)+2,2,:),1,2)]*grad_dirichlet(macMesh.node(macFemspace.bd_dof(index_x0(end)),:),idx))'*[-1; 0];
A_restricted = A(macFemspace.bd_dof(index_x0(2:end-1)),:);
residual = A_restricted*macSol;
rhs = residual - (h/6)*[flux_at_first_node; zeros(number_internal_nodes-2,1); flux_at_last_node];
flux = T\rhs;
f4 = [flux_at_first_node; flux; flux_at_last_node];

clear flux
flux = [-f1, f2, f3, -f4];