function decompose(RB) %[Theta, ThetaF, AX, Fq, Ch, ChT, S] =
%DECOMPOSE creates an affine decomposition of a problem
if exist(RB.affine.file, 'file'), return; end
if ~isfield(RB,'chol'), RB.chol = true; end

rnode = RB.rmesh.node;
elem  = RB.rmesh.elem;

%% CONSTANTS
dim = size(elem,2) - 1;
d = size(rnode,2);
R = numel(RB.G);       % number of coarse subdomains
nu0 = eye(d);           % viscosity

%% AFFINE TRANSFORMATION
[J,Ginv,nu,xi] = deal(cell(R,1));
for r=1:R
    J{r} = det(RB.G{r}); % Absolute value (in general) - but here it is ok.
    Ginv{r} = inv(RB.G{r});
    
    % coefficients of the transformed system
    nu{r} = simplify(Ginv{r} * nu0 * transpose(Ginv{r}) * J{r});
    xi{r} = simplify(Ginv{r} * J{r});
end

%% affine decomposition
% We decompose partial matrices / vectors of the problem
fprintf('Affine decomposition of partial (STI, DIV, LAG, RHS) matrices\n');

[ia, ja, ra, ib, jb, rb] = deal(zeros(0,1));
[ThetaA, ThetaB, ThetaF, ThetaL] = deal(cell(1,0));
[QA, QB] = deal(0);
QF = R;
QL = R;
for r=1:R
    for i=1:d
        for j=1:d
            if (~(nu{r}(i,j) == 0)) && (i<=j)
                QA = QA + 1;
                ia(QA) = i;
                ja(QA) = j;
                ra(QA) = r;
                ThetaA{QA} = nu{r}(i,j);
            end
            if ~(xi{r}(i,j) == 0)
                QB = QB + 1;
                ib(QB) = j;
                jb(QB) = i;
                rb(QB) = r;
                ThetaB{QB} = xi{r}(i,j); % TRANSPOSED FROM ARTICLE !!!
            end
        end
    end
    ThetaF{r} = J{r};
    ThetaL{r} = J{r};
end


%% REFERENCE MESH
% reference mesh initialization
fprintf('Reference mesh - FE spaces\n');

[rmesh, father] = RB.getfine(RB);

% division into coarse subdomains
subdomain = cell(R,1);
for r=1:R
    subdomain{r} = find(father==r);
end

% finite elmenent spaces
rmesh = gradbasis(rmesh);
rfemspace = get_femspace(rmesh, RB.elemtype);
rpfemspace = get_femspace(rmesh, RB.pelemtype);

% add fields using init_poisson
rvp = struct('elemtype', rfemspace.elemtype, 'a', 1, ...
    'f', 1, 'bc', 'zero_dirichlet');
[rvp, rmesh] = init_poisson(rmesh, rfemspace, rvp);

% constants
NA = rfemspace.ndof;
NB = rpfemspace.ndof; 
NX = dim*NA + NB + 1;

% BOUNDARY
bdidx = zeros(NA,1); % boundary conditions
allfun = @(x)(true(size(x,1),1));
allflag = get_bdflag(rmesh,allfun);
bddof = get_bddof(rmesh,rfemspace,allflag);
bdidx(bddof) = 1;
TAbd = spdiags(bdidx, 0, NA, NA);
TA = spdiags(1-bdidx, 0, NA, NA);
TX = stack_blkdiag(TA,dim);

%% PARTIAL STIFFNESS MATRICES (symmetric!)
fprintf('Assembling partial stiffness matrices\n');

Aq = cell(QA,1);
Aq(:) = {sparse(NA, NA)};

[lambda,weight] = quadpts(dim, 2*(rfemspace.deg-1)); 
LA = rfemspace.ldof;

for qa=1:QA
    for i=1:LA
        for j=1:LA
            if  (j<i) %symmetric
                continue;
            end
            Aij = zeros(numel(subdomain{ra(qa)}),1);
            if ia(qa) == ja(qa)
                loci = ia(qa);
                locj = ja(qa);
            else
                loci = [ia(qa), ja(qa)];
                locj = [ja(qa), ia(qa)];
            end
            for m=1:numel(weight)
                for t=1:numel(loci)
                    Aij = Aij + weight(m) .* ...
                        evalb(rmesh, subdomain{ra(qa)}, ...
                        lambda(m,:), i, loci(t), ...
                        rfemspace.elemtype) .* ...
                        evalb(rmesh, subdomain{ra(qa)}, ...
                        lambda(m,:), j, locj(t), ...
                        rfemspace.elemtype);
                end
            end
            Aij = Aij .* rmesh.volume(subdomain{ra(qa)});
            Aq{qa} = Aq{qa} + sparse(...
                double(get_dof(rmesh, subdomain{ra(qa)}, rfemspace, i)), ...
                double(get_dof(rmesh, subdomain{ra(qa)}, rfemspace, j)), ...
                Aij, NA, NA);
            if  (i<j)
                Aq{qa} = Aq{qa} + sparse(...
                    double(get_dof(rmesh, subdomain{ra(qa)}, ...
                    rfemspace, j)), ...
                    double(get_dof(rmesh, subdomain{ra(qa)}, ...
                    rfemspace, i)), ...
                    Aij, NA, NA);                
            end
        end
    end
    Aq{qa} = TA * Aq{qa} * TA; % Dirichlet boundary conditions
end

% add the last one - enforcing the Dirichlet BC
QA = QA + 1;
Aq{QA} = TAbd;
ThetaA{QA} = sym(1);

%% PARTIAL DIV MATRICES 
fprintf('Assembling partial div matrices\n');

Bq = cell(QB,1);
Bq(:) = {sparse(NB, dim*NA)};

[lambda, weight] = quadpts(dim, rfemspace.deg-1+rpfemspace.deg);
LB = rpfemspace.ldof;

for qb = 1:QB
    for j=1:LA
        for k=1:LB
            Bij = zeros(numel(subdomain{rb(qb)}),1);
            for i= 1:numel(weight)
                Bij = Bij + weight(i) * ...
                    evalb(rmesh, subdomain{rb(qb)}, lambda(i,:), ...
                    j, jb(qb), rfemspace.elemtype) .* ...
                    evalb(rmesh, subdomain{rb(qb)}, lambda(i,:), ...
                    k, 0,  rpfemspace.elemtype);
            end
            Bij = Bij .* rmesh.volume(subdomain{rb(qb)});
            Bq{qb} = Bq{qb} - sparse(... % MINUS !!!!
                double(get_dof(rmesh,subdomain{rb(qb)},rpfemspace,k)), ...
                double(get_dof(rmesh,subdomain{rb(qb)},rfemspace,j)) ...
                + (ib(qb)-1)*NA, ...
                Bij, NB, dim*NA);
        end
    end
    Bq{qb} = Bq{qb} * TX; % Dirichlet boundary conditions
end

%% PARTIAL LAGRANGE MULTIPLIER FIELDS
fprintf('Assembling partial multiplier fields\n');

Lq = cell(QL,1);
Lq(:) = {sparse(NB, 1)};

[lambda, weight] = quadpts(dim, rpfemspace.deg);

for ql = 1 :QL
    for j=1:LB
        lag = 0;
        for i=1:numel(weight)  
            lag = lag + weight(i) * evalb(rmesh, subdomain{ql}, ...
                lambda(i,:), j, 0, rpfemspace.elemtype);
        end
        lag = lag .* rmesh.volume(subdomain{ql});
        Lq{ql} = Lq{ql} + accumarray(get_dof(rmesh, subdomain{ql}, ...
            rpfemspace, j), lag, [NB, 1]);
    end
end

%% PARTIAL RHS VECTORS (FULL SIZE)
fprintf('Assembling partial RHS\n');

Fq = cell(QL,d);
Fq(:) = {zeros(NA, dim)};
[lambda, weight] = quadpts(dim, rfemspace.deg);

for out=1:d
    fun = str2func(['fstokesmicro' num2str(out)]);
    for qf = 1:QF
        for j=1:LA
            arhs = 0;
            for i=1:numel(weight)  % L2 PART OF FORCE FIELD
                eval_fun = evalf(rmesh, subdomain{qf}, [], lambda(i,:), fun);
                arhs = arhs + weight(i) * ...
                    bsxfun(@times, eval_fun, evalb(rmesh, subdomain{qf}, ...
                    lambda(i,:), j, 0, rfemspace.elemtype));
            end
            arhs = bsxfun(@times, arhs, rmesh.volume(subdomain{qf}));
            for i=1:dim
                Fq{qf,out}(:,i) = Fq{qf,out}(:,i) + ...
                    accumarray(get_dof(rmesh, subdomain{qf}, rfemspace, j), arhs(:,i), [NA, 1]);
            end
        end
        Fq{qf,out}(bddof,:) = 0; % Dirichlet boundary conditions
        Fq{qf,out} = reshape(Fq{qf,out},[],1);
        Fq{qf,out}(NX) = 0;
    end
end

%% AFFINE DECOMPOSITION OF THE FULL SYSTEM
fprintf('Creating affine decomposition of the full system\n');

Q = QA + QB + QL;
Theta = cell(Q,1);
AX = cell(Q,1);
for qa = 1:QA
    Theta{qa} = ThetaA{qa};
    AX{qa} = [stack_blkdiag(Aq{qa},dim), sparse(dim*NA,NB+1);
        sparse(NB+1,dim*NA), sparse(NB+1,NB+1)];
end
for qb = 1:QB
    Theta{qb+QA} = ThetaB{qb};
    AX{qb+QA} = [sparse(dim*NA,dim*NA), transpose(Bq{qb}), sparse(dim*NA,1);
        Bq{qb}, sparse(NB,NB+1);
        sparse(1,dim*NA+NB+1)];
end
for ql = 1:QL
    Theta{ql+QA+QB} = ThetaL{ql};
    AX{ql+QA+QB} = [sparse(dim*NA,dim*NA+NB+1);
        sparse(NB,dim*NA+NB), Lq{ql};
        sparse(1,dim*NA),transpose(Lq{ql}), sparse(1,1)];
end

%% REDUCTION OF THE AFFINE REPRESENTATION
fprintf('Symbolic reduction of the affine representation\n');
if isfield(RB,'reduce')
  [Theta, ThetaF, AX, Fq] = RB.reduce(Theta, ThetaF, AX, Fq); 
end

%% FULL MASS MATRIX
fprintf('Assemble the full mass matrix (with lambda modification)\n');
% velocity (1-dim) stiffness matrix
Aref = assemble_stima(rmesh, rfemspace, rvp);
Aref = TA*Aref*TA + TAbd;
% velocity mass matrix
MV = assemble_mass(rmesh, rfemspace);
MV = TA*MV*TA;
% NORM CORRECTOR
if RB.chol
  fprintf('Eigenproblem for the norm corrector\n');
  %eigoptLA = struct('tol', 1e-4,'sigma','SA','disp',true);
  %lambda   = bleigifp(Aref,MV,1,eigoptLA);
  %lambda   = lambda(1);
  lambda = eigs(Aref,MV,1,'sm',struct('issym',true));
else
  lambda=0;
end
%pressure mass matrix
MP = assemble_mass(rmesh, rpfemspace);
% X - scalar product matrix
MX = blkdiag(stack_blkdiag(Aref+lambda*MV, dim), MP, 1); % X-scalar product


if RB.chol
  fprintf('Cholesky factorization of the mass matrix\n');
  % factorization: MX  = S * Ch' * Ch  * S'
  %                MX\ = S * Ch  \ Ch' \ S'
  [Ch,pos,S] = chol(MX);
  if pos ~= 0,
    error('Mass matrix not positive definite!');
  end
end

fprintf('Saving the values\n');
if RB.chol
  save(RB.affine.file,'Theta','ThetaF','AX','Fq','MX','Ch','S');
else
  save(RB.affine.file,'Theta','ThetaF','AX','Fq','MX','-v7.3');
end
end