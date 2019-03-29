function stats = rbStokesOffline(RB, rbSample)
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

%% initialization
stats = struct('assemble', 0, 'solve', 0, 'scm', [], 'greedy', []); %('online',0,'nonline',0,'assemble',0,'solve',0,'scm', struct, greedy, cell(0));
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
NX  = size(AX{1},1);

%% OFFLINE VARIABLES

% FINE SOLUTIONS: SOL{i}(:,n) = U^{i,n}
SOL = cell(dim,1);  SOL(:) = {zeros(NX,0)};   
% SUPREMIZERS: SUP{i}(:,q,n) = T^q(U^{i,n})
SUP = cell(dim,1);  SUP(:) = {zeros(NX,Q,0)};
% F-SUPREMIZERS: EF{i}(:,q) = G^{i,q}
EF  = cell(dim,1);  EF(:)  = {zeros(NX,QF)};
% CHOSEN PARAMETERS SN, RB SIZES N
SN  = cell(dim,1);  SN(:)  = {zeros(0,P)};
N   = zeros(dim,1);

%% ONLINE VARIABLES

% OM{i}(n,m,q,r) = ( T^q( U^{i,n} ), T^r( U^{i,m} ) )_X
OM =cell(dim,1);   OM(:) = {zeros(0,0,Q,Q)};
% ON{i}(n,q,r) = F^{i,r}( T^q( U^{i,n} ) )
ON =cell(dim,1);   ON(:) = {zeros(0,Q,QF)};
% OP{i}(q,r) = F^{i,q}( G^{i,r} )
OP =cell(dim,1);   OP(:) = {zeros(QF,QF)};
% OS{i,j}(n,q) = F^{i,q}( U^{j,n} )
OS =cell(dim,dim); OS(:) = {zeros(0,QF)};
% OT{i,j}(n,m,q) = ( U^{i,n}, T^{q}( U^{j,m} ) )_X
OT =cell(dim,dim); OT(:) = {zeros(0,0,Q)};

% a posteriori coefficients for online stage - F PART
fprintf('A posteriori precomputation.\n');
assemble = tic;
if ~RB.chol
  agmg(MX,[],[],[],[],[],[],1);
end
for i=1:dim
  for qf = 1:QF
    if RB.chol
      EF{i}(:,qf) = S*(Ch\(ChT\(S'*Fq{qf,i})));
    else
      EF{i}(:,qf) = agmg(MX, Fq{qf,i},[],1e-10,[],true,[],2);
    end
  end
  OP{i} = EF{i}' * MX * EF{i};
end
stats.assemble = stats.assemble + toc(assemble);
stats.greedy = cell(dim,1);

%% GREEDY ALGORITHM FOR EACH OUTPUT
for i=1:dim
  totest = true(sampleSize,1);
  apoest = zeros(sampleSize,1);
  stats.greedy{i} = struct('maxres',[], 'online', [], 'nonline', []);
  while N(i) <= opt.Nmax
    if N(i)==0
      newi = randi(sampleSize);
    else
      apoest(totest) = rbOnline(i, totest);
      totest = apoest > opt.epstol;
      [maxres, newi] = max(apoest);
      fprintf('dir: %d/%d, BF: %d, max res: %12.15f\n',i,dim,N(i),maxres);
      stats.greedy{i}.maxres(N(i)) = maxres;
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
save(RB.online.file, 'OM','ON','OP','OS','OT','N','-v7.3');
if opt.save
  save(RB.offline.file, 'SOL','SUP','EF','SN','N','stats','-v7.3');
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
      
      apoest(tr) = ...
        TF * OP{i} * TF' - ...
        2 * (Ubar' * Fbar) + ...
        Ubar' * Abar * Ubar;
    end
    apoest(apoest<0) = 0;
    apoest = sqrt(apoest);
    if ~opt.fakescm
      apoest = apoest./betaSCM(ind);
    end
    stats.greedy{i}.nonline(N(i)) = testSize;
    stats.greedy{i}.online(N(i)) = toc(online);
  end

  function update_online_fields(i,N,newmu)
    SN{i}(N(i),1:P) = newmu; %#ok 
    
    % NEW RB SOLUTION
    fprintf('solve / ');
    solve = tic;
    if opt.refdom
      SOL{i}(:,N(i)) = assemble_affine(newmu, Theta, AX) \ ...
        assemble_affine(newmu, ThetaF, Fq(:,i));
    else
      mesh = get_transformed_domain(RB, newmu);      
      vp = struct('elemtype',RB.elemtype,'pelemtype',RB.pelemtype,...
        'a',1,'f',str2func(['fstokesmicro' num2str(i)]),'solver','iuzawapcg');
      vp.iuzawapcgopt.verbose = 2;
      vp.bc = cell(2,1);
      vp.bc{1} = [0,0,0]; vp.bc{2} = [0,0,0];
      [usol, ~, psol] = stokes(mesh, vp);
      SOL{i}(:,N(i)) = [usol(:); psol; 0];
    end
    stats.solve = stats.solve + toc(solve);
    
    % orthonormalize
    fprintf('GS / ');
    assemble = tic;
    SOL{i}(:,N(i)) = SOL{i}(:,N(i)) - SOL{i}(:,1:N(i)-1) * ...
      (SOL{i}(:,1:N(i)-1)' * MX * SOL{i}(:,N(i)));
    SOL{i}(:,N(i)) = SOL{i}(:,N(i)) / ...
      sqrt((SOL{i}(:,N(i))' * MX * SOL{i}(:,N(i))));
    
    % partial supremizers
    fprintf('supremizers. / ');
    for q=1:Q
      SUP{i}(:,q,N(i)) = AX{q}*SOL{i}(:,N(i));
    end
    if RB.chol
      SUP{i}(:,:,N(i)) = S*(Ch\(ChT\(S'*SUP{i}(:,:,N(i)))));
    else
      agmg(MX,[],[],[],[],[],[],1);
      for q=1:Q
        SUP{i}(:,q,N(i)) = agmg(MX, SUP{i}(:,q,N(i)),[],1e-10,[],true,[],2);
      end
    end
    
    % OM 
    fprintf('OM / ');
    for n=1:N(i)
      addm = SUP{i}(:,:,N(i))' * MX * SUP{i}(:,:,n);
      OM{i}(N(i),n,:,:) = reshape(addm, 1,1,Q,Q);
      OM{i}(n,N(i),:,:) = reshape(addm',1,1,Q,Q);
    end
    
    % ON
    fprintf('ON\n');
    for qf1=1:QF
      ON{i}(N(i),1:Q,qf1) = Fq{qf1,i}' * SUP{i}(:,:,N(i));
    end
    
    stats.assemble = stats.assemble + toc(assemble);
  end

  function output_eval_fields
    fprintf('OS and OT\n');
    assemble = tic;
    for ii=1:dim
      for jj=1:dim        
        OS{ii,jj} = SOL{ii}' * MX * EF{jj}; %#ok 
        for nj = 1:N(jj)
          OT{ii,jj}(1:N(ii),nj,1:Q) = SOL{ii}' * MX * SUP{jj}(:,:,nj); %#ok          
        end
      end
    end
    stats.assemble = stats.assemble + toc(assemble);
  end
end