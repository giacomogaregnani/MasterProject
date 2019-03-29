function [ustokes,ufemspace,pstokes,pfemspace,phidarcy,phifemspace,mesh,vp] = darcystokesv3(mesh, interface, vp)

%% 1. Stokes part

%% 1.1 FINITE ELEMENTS INITIALIZATION
mesh{1} = gradbasis(mesh{1});
femspace{1} = get_femspace(mesh{1}, vp{1}.elemtype);
pfemspace = get_femspace(mesh{1}, vp{1}.pelemtype);
[vp{1}, mesh{1}] = init_stokes(mesh{1}, femspace{1}, vp{1});

%% 1.2 CONSTANTS
dim = size(mesh{1}.elem, 2)-1;
NA = femspace{1}.ndof;
NB = pfemspace.ndof; 

%% 1.3 ASSEMBLE STIFFNESS MATRIX
AS = assemble_stima(mesh{1}, femspace{1}, vp{1});
AS = stack_blkdiag(AS,dim);

%% 1.4 ASSEMBLE RIGHT HAND SIDE
fS = assemble_rhs(mesh{1}, femspace{1}, vp{1});
fS = reshape(fS, [dim*NA, 1, vp{1}.NF]);

%% 1.5 ASSEMBLE DIV MATRIX
BS = assemble_div(mesh{1}, femspace{1}, pfemspace, vp{1}); 
gS = zeros(NB, 1, vp{1}.NF);

%% 2. Darcy part

%% 2.1 FE INITIALIZATION
mesh{2} = gradbasis(mesh{2});
femspace{2} = get_femspace(mesh{2}, vp{2}.elemtype);
[vp{2}, mesh{2}] = init_poisson(mesh{2}, femspace{2}, vp{2});

%% 2.2 ASSEMBLE STIFFNESS MATRIX 
AD = assemble_stima(mesh{2}, femspace{2}, vp{2});
NAD = size(AD,1);

%% 2.2 ASSEMBLE RHS 
fD = assemble_rhs(mesh{2}, femspace{2}, vp{2});


%% 3. INTERFACE COMPONENTS

%% 3.1 symmetric condition
%ASYM = assemble_adsym;
ASYM = assemble_interface_sym(mesh,femspace,vp);
%% 3.2 non-symmetric condition
ANSYM = assemble_interface_nonsym(mesh,femspace,interface);

%% assemble the whole system
A11 = vp{2}.n .* AS + vp{2}.n .* vp{2}.alpha.*ASYM; %for 3d problems, alpha could be direction-dependent.
A12 = vp{2}.n .* vp{2}.g .* vp{2}.rho.*ANSYM;
A22 = vp{2}.g .* vp{2}.rho .* AD;
A31 = vp{2}.n .* BS;

f1 = vp{2}.n.*fS;
f2 = vp{2}.g.*vp{2}.rho.*fD;
f3 = vp{2}.n.*gS;

%% boundary conditions (only zero Dirichlet or zero Neumann)
boundary_conditions;

A = [A11, A12, A31'; 
  -A12', A22, sparse(NAD,NB); 
  A31, sparse(NB,NAD), sparse(NB,NB)]; 
rhs = [f1; f2; f3];
 
%TODO: more general formulation. Does the last term come with an vp{2}.n
%coefficient? It is now zero for all testing purposes.

%rhs =[vp{2}.n.*fS; 10.*fD; vp{2}.n.*gS]; 
% rhs =[vp{2}.n.*fS; vp{2}.g.*vp{2}.rho.*fD; vp{2}.n.*gS]; 

%% solve
%Use backslash solver for the 2d system
Sols = A\rhs;
%Split the solution vector in the three relevant components ustokes,
%phidarcy, pstokes

ustokes = reshape(Sols(1:dim*NA,:), NA, dim, []); %recover a vector field
phidarcy = Sols((size(AS,1)+1):(size(AS,1)+size(AD,1))); % Darcy discharge
pstokes = Sols((size(AS,1)+size(AD,1)+1):end); % Stokes pressure

%Do we want a 0 mean pressure? Would make sense if we have Dirichlet bc on
%the exterior
%pstokes = bsxfun(@minus, pstokes, get_average(mesh{1}, pfemspace, pstokes));
        
%% Define femspaces to be returned by the function
ufemspace = femspace{1};
phifemspace = femspace{2};

  function boundary_conditions
    %% STOKES part 1
    % FIND DIRICHLET DOF
    vp{1}.ispureneumann   = numel(mesh{1}.bdflag{1}) == 0;
    vp{1}.ispuredirichlet = numel(mesh{1}.bdflag{2}) == 0;
    if vp{1}.ispureneumann
      bddof1 = get_dof(mesh{1}, 1, femspace{1}, 1);
      bdcoor1 = get_rc(mesh{1}, [], 1, [], femspace{1}.nodelambda(1,:));
    else
      [bddof1, bdcoor1] = get_bddof(mesh{1}, femspace{1}, mesh{1}.bdflag{1});
    end
    
    % boundary "logical" field
    bdidx1 = zeros(dim*NA,1);
    for j=0:dim-1
      bdidx1(NA*j + bddof1) = 1;
    end
    
    % insert Dirichlet values
    uD1 = zeros(dim*NA,1);
    if ~vp{1}.ispureneumann
      if isnumeric(vp{1}.bc{1})
        for j=1:dim
          uD1((j-1)*NA + bddof1) = vp{1}.bc{1}(j);
        end
      else
        dirval = vp{1}.bc{1}(bdcoor1);
        for j=1:dim
          uD1((j-1)*NA + bddof1) = dirval(:,j);
        end
      end
    end    
    
    % UPDATE RHS BASED ON THE DIRICHLET BC
    f1 = bsxfun(@minus,f1,A11*uD1);
    f2 = bsxfun(@minus,f2,-A12'*uD1);
    f3 = bsxfun(@minus,f3,A31*uD1);    
    
    Tbd1 = spdiags(bdidx1, 0, dim*NA, dim*NA);
    T1 = spdiags(1-bdidx1, 0, dim*NA, dim*NA);
    
    %% DARCY part 1    
    vp{2}.ispureneumann = numel(mesh{2}.bdflag{1}) == 0;
%    if vp{2}.ispureneumann
%      bddof2 = get_dof(mesh{2}, 1, femspace{2}, 1);
%      bdcoor2 = get_rc(mesh{2}, [], 1, [], femspace{2}.nodelambda(1,:));
%    else
      [bddof2, bdcoor2] = get_bddof(mesh{2},femspace{2},mesh{2}.bdflag{1});
%    end
    
    % boundary "logical" field
    bdidx2 = zeros(NAD,1);
    bdidx2(bddof2) = 1;
    
    % (ALMOST EVERYWHERE) ZERO SOLUTION WITH RIGHT DIRICHLET BC
    uD2 = zeros(NAD,1);
    if ~vp{2}.ispureneumann
      if isnumeric(vp{2}.bc{1})
        if vp{2}.bc{1} ~= 0
          uD2(bddof2) = vp{2}.bc{1};
        end
      else
        uD2(bddof2) = vp{2}.bc{1}(bdcoor2);
      end
    end
    
    % UPDATE RHS BASED ON THE DIRICHLET BC
    f2 = bsxfun(@minus, f2, AD*uD2);
    
    Tbd2 = spdiags(bdidx2, 0, NAD, NAD);
    T2 = spdiags(1-bdidx2, 0, NAD, NAD);
    
    %% update matrices
    
    A22 = T2*A22*T2+ Tbd2;
    A11 = T1*A11*T1+ Tbd1;
    A31 = A31*T1;
    A12 = T1*A12*T2;
    
    % ASSUME NEUMANN B.C. are zero
    
    % ENFORCE STRONGLY DIRICHLET B.C.
    f1(bdidx1==1) = uD1(bdidx1==1);
    f2(bdidx2==1) = uD2(bdidx2==1);
    
  end

end