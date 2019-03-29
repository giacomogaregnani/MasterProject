function [usol,ufemspace,psol,pfemspace,phisol,phifemspace,mesh,vp] = darcystokesv2(mesh, interface, vp, n, alpha)

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
ASB = assemble_adsym;

%% 3.2 non-symmetric condition
ANS = assemble_adnonsym;



%% assemble the whole system
%Check that dimensions are consistent
% sA111 = size(AS)
% sA112 = size(ASB)
% sA12 = size(ANS)
% sA13 = size(BS')
% sA21 = size(ANS')
% sA22 = size(AD)
% sA31 = size(BS)
% sfs = size(fS)
% sfd = size(fD)

A11 = n*AS + n*alpha*ASB;
A12 = n*ANS;
A22 = AD;
A31 = n*BS;

%% boundary conditions (only zero Dirichlet or zero Neumann)
boundary_conditions;

A = [A11, A12, A31'; 
  -A12', A22, sparse(NAD,NB); 
  A31, sparse(NB,NAD), sparse(NB,NB)]; 
% sa = size(A)
%rhs = [n.*fS; fD; zeros(size(BS,1),1)]; 
%TODO: more general formulation. Does the last term come with an n
%coefficient? It is now zero for all testing purposes.
rhs =[fS; fD; gS]; 

%% solve
%Use backslash solver for the 2d system
Sols = A\rhs;
%Split the solution vector in the three relevant components u_f, phi_p, p_f
%usol = Sols(1:size(AS,1));
%phisol = Sols((size(AS,1)+1):(size(AS,1)+size(AD,1)));
%psol = Sols((size(AS,1)+size(AD,1)+1):end);

usol = reshape(Sols(1:dim*NA,:), NA, dim, []);
phisol = Sols((size(AS,1)+1):(size(AS,1)+size(AD,1)));
psol = Sols((size(AS,1)+size(AD,1)+1):end);

%Do we want a 0 mean pressure? Would make sense if we have Dirichlet bc on
%the exterior
psol = bsxfun(@minus, psol, get_average(mesh{1}, pfemspace, psol));
        

%%Define femspaces to be returned by the function
ufemspace = femspace{1};
phifemspace = femspace{2};

%% disassemble the whole system to separate functions in function spaces


  function ASB = assemble_adsym
    ASB = sparse(dim*NA, dim*NA);
    % quadrature formula on edges
    [lambda, weight] = quadpts(dim-1, 2*femspace{1}.deg);
    NQ = numel(weight);
    % for each js = 1:dim we consider those interface elements that have the
    % boundary edge oppoiste to the node with local coordinate js. This
    % simplifies how we obtain the degrees of freedom.
    for js=1:dim+1
      ind = mesh{1}.bdflag{3}(:,2) == js;
      % in this pass we only work with interface elements with numbers swhe
      swhe = mesh{1}.bdflag{3}(ind, 1);
      if isempty(swhe),
        continue;
      end
      % compute normals and tangents at the considered elements. Here we assume
      % that the dimension is 2. One of the few parts that will need
      % generalization before coming to 3D
      normals = get_normals(mesh{1}, swhe, js);
      tangents = [-normals(:,2), normals(:,1)]; % works only in 2D
      % modify the quadrature formula (just the barycentric coordinates) to
      % simplify the integration: we integrate over elements.
      lam = [lambda(:,1:js-1), zeros(size(lambda,1),1), lambda(:,js:dim)];
      % which local degrees of freedom are of interest (those that lie on the
      % interface edge).
      ldof = find(femspace{1}.nodelambda(:,js)==0)';
      % volumes (lengths in 2D, surfaces in 3D) of the interface edges.
      ssvol = subsimplex_volume(mesh{1}, [], swhe, [1:js-1,js+1:dim+1]);
      for k=1:NQ
        % coef is the part alpha * n / sqrt(tau * K * tau), computed at the
        % quadrature points. We assume that we are in 2 dimensions so that
        % there is just one tau.
        coef =0;
        % if the tensor is a function handle, we will need coordinates of the
        % quadrature nodes
        if ~isnumeric(vp{2}.a)
          xloc = get_rc(mesh{1},[],swhe,[],lam(k,:));
        end
        % now we assemble the values tau_{t1} * K_{t1, t2} * rau_{t2} and sum
        % them up
        for t1 = 1:dim
          for t2 = 1:dim
            if isnumeric(vp{2}.a)
              if numel(vp{2}.a) == 1
                if t1==t2,
                  vpa = vp{2}.a;
                else
                  continue;
                end
              else
                vpa = vp{2}.a(t1, t2);
                if vpa == 0
                  continue;
                end
              end
            else
              vpa = vp{2}.a(xloc,t1,t2);
            end
            coef = coef + vpa .* tangents(:,t1) .* tangents(:,t2);
          end
        end
        % the assembled coef is now modified to its final form, where
        % alpha is the Beavers-Saffman constant determined experimentally.
        coef = 1./sqrt(coef);
        % we continue by assembling the required matrix. We loop over local
        % degrees of freedom (of interest) twice and then assemble the values
        % into a sparse matrix.
        for m1=ldof
          % degrees of freedom (dof1) and basis functions (evb1) evaluated at
          % the quadrature points
          dof1 = get_dof(mesh{1}, swhe, femspace{1}, m1);
          evb1 = evalb(mesh{1}, swhe, lam(k,:), m1, 0, femspace{1}.elemtype);
          for m2=ldof
            % degrees of freedom (dof2) and basis functions (evb2) evaluated at
            % the quadrature points
            dof2 = get_dof(mesh{1}, swhe, femspace{1}, m2);
            evb2 = evalb(mesh{1}, swhe, lam(k,:), m2, 0, femspace{1}.elemtype);
            for t1 = 1:dim
              for t2 = 1:dim
                Aij = evb1 .* tangents(:,t1) .* evb2 .* tangents(:,t2) .* ...
                  weight(k) .* ssvol .* coef;
                ASB = ASB + sparse(double(dof1+(t1-1)*NA), double(dof2+(t2-1)*NA), Aij, dim*NA, dim*NA);
              end
            end
          end
        end
      end
    end
  end

  function BS = assemble_adnonsym
    BS = sparse(dim*NA, NAD);
    % quadrature formula on edges
    [lambda, weight] = quadpts(dim-1, femspace{1}.deg + femspace{2}.deg);
    NQ = numel(weight);
    dof = cell(2,1);
    % We go through interface{1} and interface{2}, one at a time.
    for r = 1:2
      % for each js = 1:dim we consider those interface elements in mesh{r} 
      % that have the boundary edge oppoiste to the node with local 
      % coordinate js. This simplifies how we obtain the degrees of freedom
      for js=1:dim+1
        ind = mesh{r}.bdflag{3}(:,2) == js;
        % in this pass we only work with interface elements with numbers
        % whe1 (in mesh{r})
        whe1 = mesh{r}.bdflag{3}(ind, 1);
        % we restrict to those elements that are present in 
        % interface{r}(:,1)
        [~,i1,i2] = intersect(whe1, interface{r}(:,1));
        whe1 = whe1(i1);
        % corresponding interface elements from mesh{3-r} have element
        % numbers whe2
        whe2 = interface{r}(i2,2);
        if isempty(whe1), continue; end
        % compute normals 
        normals = (-1)^(r+1)*get_normals(mesh{r}, whe1, js); % TODO: watch out for direction!
        %Antoine: 27/3, changed sign
        % modify the quadrature formula (just the barycentric coordinates) to
        % simplify the integration: we integrate over elements.
        lam1 = [lambda(:,1:js-1), zeros(size(lambda,1),1), lambda(:,js:dim)];
        % which local degrees of freedom are of interest (those that lie on the
        % interface edge).
        ldof1 = find(femspace{r}.nodelambda(:,js)==0)';
        ldof2 = 1:femspace{3-r}.ldof;
        % volumes (lengths in 2D, surfaces in 3D) of the interface edges.
        ssvol = subsimplex_volume(mesh{r}, [], whe1, [1:js-1,js+1:dim+1]);
        for k=1:NQ
          % cartesian coordinates of quadrature nodes
          xloc = get_rc(mesh{r},[],whe1,[],lam1(k,:));
          % barycentric coordinates for integration in mesh{3-r}
          lam2 = get_lc2(mesh{3-r}, whe2, xloc);          
          % We loop over local
          % degrees of freedom (of interest) in mesh{r}
          for m1=ldof1
            % degrees of freedom (dof1) and basis functions (evb1) evaluated at
            % the quadrature points in mesh{r}
            dof{r} = get_dof(mesh{r}, whe1, femspace{r}, m1);
            evb1 = evalb(mesh{r}, whe1, lam1(k,:), m1, 0, femspace{r}.elemtype);
            % We loop over local
            % degrees of freedom (of interest) in mesh{3-r}
            for m2=ldof2
              % degrees of freedom (dof2) and basis functions (evb2) evaluated at
              % the quadrature points in mesh{3-r}
              dof{3-r} = get_dof(mesh{3-r}, whe2, femspace{3-r}, m2);
              evb2 = evalb(mesh{3-r}, whe2, lam2, m2, 0, femspace{3-r}.elemtype);
              for t = 1:dim
                Aij = evb1 .* evb2 .* normals(:,t) .* weight(k) .* ssvol;
                BS = BS + sparse(double(dof{1}+(t-1)*NA), double(dof{2}), Aij, dim*NA, NAD);
              end
            end
          end
        end
      end
    end
  end


  function boundary_conditions
    %% STOKES part 1
    % FIND DIRICHLET DOF
    vp{1}.ispureneumann = numel(mesh{1}.bdflag{1}) == 0;
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
    fS = bsxfun(@minus,fS,A11*uD1);
    fD = bsxfun(@minus,fD,-A12'*uD1);
    gS = bsxfun(@minus,gS,A31*uD1);
    
    
    Tbd1 = spdiags(bdidx1, 0, dim*NA, dim*NA);
    T1 = spdiags(1-bdidx1, 0, dim*NA, dim*NA);
    
    %% DARCY part 1
    
    vp{2}.ispureneumann = numel(mesh{2}.bdflag{1}) == 0;
    if vp{2}.ispureneumann
      bddof2 = get_dof(mesh{2}, 1, femspace{2}, 1);
      bdcoor2 = get_rc(mesh{2}, [], 1, [], femspace{2}.nodelambda(1,:));
    else
      [bddof2, bdcoor2] = get_bddof(mesh{2},femspace{2},mesh{2}.bdflag{1});
    end
    
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
    fD = bsxfun(@minus, fD, AD*uD2);
    
    Tbd2 = spdiags(bdidx2, 0, NAD, NAD);
    T2 = spdiags(1-bdidx2, 0, NAD, NAD);
    
    %% update matrices
    
    
    
    
    A22 = T2*A22*T2+ Tbd2;
    A11 = T1*A11*T1+ Tbd1;
    A31 = A31*T1;
    A12 = T1*A12*T2;
    
    % ASSUME NEUMANN B.C. are zero
    
    % ENFORCE STRONGLY DIRICHLET B.C.
    fS(bdidx1==1) = uD1(bdidx1==1);
    fD(bdidx2==1) = uD2(bdidx2==1);
    
  end

end