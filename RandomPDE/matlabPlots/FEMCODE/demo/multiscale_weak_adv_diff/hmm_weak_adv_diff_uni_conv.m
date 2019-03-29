%% USER INPUT

micElem = 'p1';
Nmic_max = 2;
macElem  = 'p1';
Nmac_max = 5;
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

expName  = 'WeakAdvDiffMS2D';
dim = 2;

h1 = zeros(Nmic_max,Nmac_max);
l2 = zeros(Nmic_max,Nmac_max);
h1s = zeros(Nmic_max,Nmac_max);
dof_mac = zeros(Nmic_max,Nmac_max);

% specify DG parameters
if length(macElem) == 3 && strcmp(macElem(3),'d')
    if ~exist('diffpart','var'),   diffpart = 'iipg'; ...
        warning('IIPG for diffusion, as not specified'); end
    if ~exist('alpha','var'),   alpha = 10; ...
        warning('default penalization parameter, as not specified'); end
    dg = zeros(Nmic_max,Nmac_max);
    dgs = zeros(Nmic_max,Nmac_max);
end

%% curves for different micro discretizations
for imic = 1:Nmic_max
    
    clear micTensor micVp micMesh
    
    %% micro problem definition
    micTensor = str2func(['aeps_' expName]);
    micVector = str2func(['beps_' expName]);
    if strcmp(coupling,'periodic')
        micVp = struct('elemtype', micElem, 'a', micTensor, 'b', micVector);
        micMesh = structured_mesh(ratio*repmat([0,1],[1,dim]),2*2^imic,struct('periodic',true, 'centre', true));
    elseif strcmp(coupling,'dirichlet');
        micVp = struct('elemtype', micElem, 'bc', 0, 'a', micTensor, 'b', micVector);
        micMesh = structured_mesh(ratio*repmat([0,1],[1,dim]),2*2^imic,struct('centre',true));
    else
        error('unknown multiscale coupling');
    end
    micVp.f = permute([zeros(dim,1), -eye(dim)],[3,2,1]);
    micMesh.bdflag = 'dirichlet'; % applies to all the boundary = zero boundary

    
    %% refine macro discretization
    for imac = 1:Nmac_max
        
        clear macRhs macVp macMesh macSol macFemspace;
        
        % macro problem definition
        macRhs = str2func(['rhs_' expName]);
        macVp = struct('elemtype', macElem, 'bc', 0, 'f', macRhs);
        
        % specify DG parameters (if the case)
        if length(macElem) == 3 && strcmp(macElem(3),'d')
            macVp.diffpart = diffpart;
            macVp.alpha = alpha;
        end
        
        macMesh = structured_mesh(repmat([0,1],[1,dim]),2^imac,struct('centre',true));
        macMesh.bdflag = 'dirichlet';
        
        
        %% run HMM
        options = struct();
        
        [macSol, macFemspace, macMesh, macVp] = ...
          hmm_weak_elliptic_uni(macMesh, macVp, micMesh, micVp, options);
      
        
        % calculate error
        dof_mac(imic,imac) = macFemspace.ndof;
        if length(macElem) == 2
            [h1(imic,imac), l2(imic,imac), h1s(imic,imac)] = ...
                get_H1error_exact(macMesh, macFemspace, macSol, @ex_PoissonMS2D);
        elseif length(macElem) == 3 && strcmp(macElem(3),'d')
            [h1(imic,imac), l2(imic,imac), h1s(imic,imac), dg(imic,imac), dgs(imic,imac)] = ...
                get_DGerror_diffusion_exact(macMesh, macFemspace, macVp, macSol, @ex_PoissonMS2D);
        end
      
    end
    
end

%% plot error curves

order = str2double(macElem(2));

figure;
for imic=1:Nmic_max 
    loglog(dof_mac(imic,:),l2(imic,:),'r-x', ...
        [10,10^3],(([10,10^3]).^(-1/2)).^(order+1),'r--');
    hold on;
end
title('Error in L2-norm');
xlabel('DOF');
ylabel('absolute error');


figure;
for imic=1:Nmic_max 
    loglog(dof_mac(imic,:),h1(imic,:),'b-x', ...
        [10,10^3],5*(([10,10^3]).^(-1/2)).^order,'b--');
    hold on;
end
title('Error in H1-norm');
xlabel('DOF');
ylabel('absolute error');

if length(macElem) == 3 && strcmp(macElem(3),'d')
    figure;
    for imic=1:Nmic_max 
        loglog(dof_mac(imic,:),dg(imic,:),'g-x', ...
            [10,10^3],5*(([10,10^3]).^(-1/2)).^order,'g--');
        hold on;
    end
    title('Error in DG-norm');
    xlabel('DOF');
    ylabel('absolute error');
end
