function stats = rbStokesOffline_files(RB, rbSample)
%RB_OFFLINE Summary of this function goes here
%   Detailed explanation goes here
%
% stats is a structure with subfields:
%  online: single-threaded time spent on online computations
%  nonline: number of online calculations
%  assemble: time spent on assembling values for online output
%  solve: time spent on solving the offline problems
%  scm: structure that describes the time cost of the SCM part

opt = RB.offline;
if ~isfield(opt,'Nmax'),    opt.Nmax = 150; end
if ~isfield(opt,'epstol'),  opt.epstol = 1e-6; end
if ~isfield(opt,'fakescm'), opt.fakescm = false; end
if ~isfield(opt,'save'),    opt.save = false; end
if ~isfield(opt,'refdom'),  opt.refdom = true; end

if ~exist(opt.folder,'dir'),
  mkdir(opt.folder);
end
%% initialization
stats = struct('online',0,'nonline',0,'assemble',0,'solve',0,'scm', struct);
P     = numel(RB.param.ref);
sampleSize = size(rbSample,1);

%% LOAD AFFINE DECOMPOSITION (coefficients)
[Theta, ThetaF] = deal([]);  
load(RB.affine.file, 'Theta', 'ThetaF');
Q   = numel(Theta);
QF  = numel(ThetaF);

%% pre-evaluation of Theta fields
ThetaField  = evalTheta(rbSample, Theta);
ThetaFieldF = evalTheta(rbSample, ThetaF);

%% SCM FOR INF-SUP LOWER BOUNDS
if ~opt.fakescm
    [betaSCM, stats.scm] = natural_SCM(RB, rbSample, false);
    if any(betaSCM <= 0),
      error('SCM did not work properly');
    end
end

%% LOAD AFFINE DECOMPOSITION (terms)
[MX, Ch, S, AX, Fq] = deal([]);  %make them global
if RB.chol
  load(RB.affine.file, 'MX', 'Ch', 'S', 'AX', 'Fq');
  ChT = Ch';
else
  load(RB.affine.file, 'MX', 'AX', 'Fq');
end
dim = size(Fq,2);
%NX  = size(AX{1},1);

%% OFFLINE VARIABLES
RED = cell(dim,1);  RED(:) = {cell(0,0)};
SUP = cell(dim,1);  SUP(:) = {cell(Q,0)};
EF  = cell(dim,1);  EF(:)  = {cell(QF,0)};

SN  = cell(dim,1);  SN(:)  = {zeros(0,P)};
N   = zeros(dim,1);

%% ONLINE VARIABLES
OM =cell(dim,1);   OM(:) = {zeros(0,0,Q,Q)};
ON =cell(dim,1);   ON(:) = {zeros(0,Q,QF)};
OP =cell(dim,1);   OP(:) = {zeros(QF,QF)};
OQ =cell(dim,1);   OQ(:) = {zeros(0,Q,QF)};
OS =cell(dim,dim); OS(:) = {zeros(0,QF)};
OT =cell(dim,dim); OT(:) = {zeros(0,0,Q)};
OU =cell(dim,dim); OU(:) = {zeros(QF,0)}; % output evaluation

% a posteriori coefficients for online stage - F PART
fprintf('A posteriori precomputation.\n');
assemble = tic;
if ~RB.chol
  agmg(MX,[],[],[],[],[],[],1);
end
for i=1:dim
  for qf = 1:QF
    EF{i}{qf} = gen_file('EF',i,qf);
    if RB.chol
      EF{i}{qf}.m = S*(Ch\(ChT\(S'*Fq{qf,i})));
    else
      EF{i}{qf}.m = agmg(MX, Fq{qf,i}, [], 1e-10, [], true, [], 2);
    end
  end
  for qf = 1:QF
    llef = EF{i}{qf}.m;
    for qff = 1:QF
      OP{i}(qf,qff) = llef' * MX * EF{i}{qff}.m;
    end
  end
end
stats.assemble = stats.assemble + toc(assemble);


%% GREEDY ALGORITHM FOR EACH OUTPUT
for i=1:dim
  totest = true(sampleSize,1);
  apoest = zeros(sampleSize,1);
  while N(i) <= opt.Nmax
    if N(i)==0
      newi = randi(sampleSize);
    else
      apoest(totest) = rbOnline(i, totest);
      totest = apoest > opt.epstol;
      [maxres, newi] = max(apoest);
      fprintf('dir: %d/%d, BF: %d, max res: %12.15f\n',i,dim,N(i),maxres);
      if ~any(totest)
        fprintf('Desired precision reached.\n');
        break;
      end
    end
    N(i) = N(i) + 1;
    newmu = rbSample(newi,:);
    update_online_fields(i,N,newmu);    
  end
  %% assign values to tensor and to errF
end
output_eval_fields;

%% save offline and online
save(RB.online.file, 'OM','ON','OP','OQ','OS','OT','OU','N','-v7.3');
if opt.save
  save(RB.offline.file, 'SN','N','stats','-v7.3');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  function apoest = rbOnline(i, totest)
    online = tic;
    ind = find(totest);
    testSize = numel(ind);
    apoest = zeros(testSize,1);
    fprintf('Solving %d online micro problems\n',testSize);
    for tr = 1:testSize
      T  = ThetaField(ind(tr),:);
      TF = ThetaFieldF(ind(tr),:);
      Abar = reshape(reshape(reshape(OM{i},[],Q)  * T', [],Q) * T', [N(i),N(i)]);
      Fbar = reshape(reshape(reshape(ON{i},[],QF) * TF',[],Q) * T', [N(i),1]);
     
      Ubar = Abar \ Fbar;
      
      apoest(tr) = TF * OP{i} * TF' + ...
        2*((Ubar' * reshape(reshape(OQ{i},[],QF) * TF',N(i),Q)) * T') + ...
        Ubar' * reshape(reshape(reshape(OM{i},[],Q) * T', [], Q) * T', N(i), N(i)) * Ubar;
    end
    apoest(apoest<0) = 0;
    apoest = sqrt(apoest);
    if ~opt.fakescm
      apoest = apoest./betaSCM(ind);
    end
    stats.nonline = stats.nonline + sampleSize;
    stats.online = stats.online + toc(online);
  end

  function update_online_fields(i,N,newmu)
    SN{i}(N(i),1:P) = newmu;
    
    % new RB function
    fprintf('Full sys. / ');
    solve = tic;
    if opt.refdom
      sol = assemble_affine(newmu, Theta, AX) \ ...
        assemble_affine(newmu, ThetaF, Fq(:,i));
    else
      mesh = get_transformed_domain(RB, newmu);
      
      vp = struct('elemtype',RB.elemtype,'pelemtype',RB.pelemtype,...
        'a',1,'f',str2func(['fstokesmicro' num2str(i)]),'solver','iuzawapcg');
      vp.iuzawapcgopt.verbose = 2;
      vp.bc = cell(2,1);
      vp.bc{1} = [0,0,0]; vp.bc{2} = [0,0,0];
      [usol, ~, psol] = stokes(mesh, vp);
      sol = [usol(:); psol; 0];
    end
    stats.solve = stats.solve + toc(solve);
    
    % orthonormalize
    fprintf('Orthonorm. / ');
    assemble = tic;
    for j=1:N(i)-1
      sol = sol - RED{i}{j}.m * (RED{i}{j}.m' * MX * sol);
    end
    sol = sol / sqrt(sol' * MX * sol);
    RED{i}{N(i)} = gen_file('RED',i,N(i));
    RED{i}{N(i)}.m = sol;
    
    % partial supremizers
    fprintf('Suprem. / ');
    if ~RB.chol
      agmg(MX,[],[],[],[],[],[],1);
    end
    for q=1:Q
      SUP{i}{q,N(i)} = gen_file('SUP',i,q,N(i));
      if ~RB.chol
        SUP{i}{q,N(i)}.m = agmg(MX, AX{q}*sol,[],1e-10,[],false,[],2);
      else
        SUP{i}{q,N(i)}.m = S*(Ch\(ChT\(S'* AX{q}*sol )));
      end
    end

    % OM - stiffness coefficients for online stage
    fprintf('Stiff. / ');
    for n=1:N(i)
      addm = zeros(Q,Q);
      lsol = RED{i}{n}.m;
      for q1 = 1:Q
        lsup = SUP{i}{q1,N(i)}.m;
        for q2 = 1:Q
          addm(q1,q2) = lsup' * (AX{q2} * lsol);
        end
      end
      OM{i}(N(i),n,:,:) = reshape(addm, 1,1,Q,Q);
      OM{i}(n,N(i),:,:) = reshape(addm',1,1,Q,Q);
    end
    
    % ON - rhs coefficients for online stage
    fprintf('RHS / ');
    for qf1=1:QF
      for q=1:Q
        ON{i}(N(i),q,qf1) = Fq{qf1,i}' * SUP{i}{q,N(i)}.m;
      end
    end
    
    % OQ, OR, a posteriori coefficients for online stage
    for qf1 = 1:QF  
      lef = EF{i}{qf1}.m;
      for q=1:Q
        OQ{i}(N(i),q,qf1) = - lef' * (AX{q} * sol);
      end
    end
    stats.assemble = stats.assemble + toc(assemble);
  end

  function output_eval_fields
    fprintf('Assembling output evaluation fields.\n');
    assemble = tic;
    for ii=1:dim
      for jj=1:dim
        for  n=1:N(ii)
          sol = RED{ii}{n}.m;
          for qf1 = 1:QF
            OS{ii,jj}(n,qf1) = sol' * MX * EF{jj}{qf1}.m;
            OU{jj,ii}(qf1,n) = Fq{qf1,jj}' * sol;
          end
          for nj=1:N(jj)
            lsol = RED{jj}{nj}.m;
            for q=1:QF
              OT{ii,jj}(n,nj,q) = - sol' * AX{q} * lsol;
            end
          end
        end    
      end
    end
    stats.assemble = stats.assemble + toc(assemble);
  end

  function m=gen_file(varargin)
    fname = RB.offline.folder;
    for ii=1:nargin
      if ischar(varargin{ii})
        fname = [fname, varargin{ii}];
      else
        fname = [fname, num2str(varargin{ii})];
      end
      if ii < nargin
        fname = [fname, '_'];
      else
        fname = [fname, '.mat'];
      end
    end
    m = matfile(fname,'writable',true);
  end
end