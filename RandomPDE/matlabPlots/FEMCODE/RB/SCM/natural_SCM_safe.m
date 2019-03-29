function [betaSCM, stats] = natural_SCM_safe(RB, M)
%NATURAL_SCM SUCCESSIVE CONSTRAINT METHOD IN NATURAL NORM with safety
%
% [betaSCM, stats] = natural_SCM(RB, scmTrain, train)
% performs the successive constraint method in natural norm with a reduced
% basis setting described in the structure RB, on a training parameter set
% given in scmTrain. If the third input parameter is set to true, we
% perform the training. If not, we just evaluate the lower bounds of the
% inf-sup constants in betaSCM. The field stats returns time costs of
% different parts:
%
% stats.eig = array of running times (in seconds) of single threaded
% evaluations of all the eigenproblems.
%
% stats.linprog = cummulative time of (in seconds) of all the solved linear
% programs.
%
% stats.nlinprog = number of solved linear programs.
%
% issues:
%   only works in 2D now.
error('This is not working, I just started implementing it but it seems to be more difficult than I thought.')

if numel(RB.param.min) ~= 2
  error('Unsupported number of parameters');
end



opt = RB.scm;
if ~isfield(opt,'theta'),  opt.theta = 0.5;   end
if ~isfield(opt,'epstol'), opt.epstol = 0.5;  end
if ~isfield(opt,'eigtol'), opt.eigtol = 1e-4; end
theta = opt.theta;

sampleSize = size(scmSample,1);
tryind = true(sampleSize,1);
Theta = [];
load(RB.affine.file, 'Theta');
Q = numel(Theta);
[betaLB, betaUB] = deal(zeros(sampleSize,0)); % initialization
ThetaVal = evalTheta(scmSample, Theta);       % PRE-evaluation of Theta Fields

%% LOAD + INIT + SETTINGS
[Ch, S, AX] = deal([]);
load(RB.affine.file, 'Ch','S','AX');
ChT = Ch';
NX = size(AX{1},1);
eigoptSA = struct('tol', opt.eigtol, 'size',NX,'sigma','SA');
eigoptLA = struct('tol', opt.eigtol, 'size',NX,'sigma','LA','disp',0);

%% INITIALIZATION (gamma, first eigenvalue)
[zvec, beta, gamma, beta2, LPA, mubar, muhat] = deal([]);
stats = struct('eig',zeros(0,1),'linprog',0,'nlinprog',0);
if ~exist(RB.scm.file,'file')
  stats = struct('eig',zeros(0,1),'linprog',0,'nlinprog',0);
  gamma = get_gamma;  k=[];
  indA1 = randi(sampleSize); minmaxUB = -Inf; % first mubar and muhat - random  
else
  load(RB.scm.file, 'zvec', 'beta', 'gamma', 'beta2', 'LPA', 'mubar', 'muhat');
  load(RB.scm.file, 'stats');
  k = 1:numel(beta);
end

%% SCM MAIN ITERATION
for zz=1:sampleSize % this will end much sooner
  if ~isempty(k)
    updateUBLB(zvec, beta, gamma, beta2, LPA, k);

    tryind = max(betaLB / theta, [], 2) <= opt.epstol;
    [~, indA1] = min(max(betaUB / theta, [], 2)); %smallest betaUB
    [~, indA2] = max(betaUB(indA1,:) / theta);
    minmaxUB   = betaUB(indA1, indA2) / theta;
    [~, indB1] = min(max(betaLB / theta, [], 2)); %smallest betaLB
    [~, indB2] = max(betaLB(indB1,:) / theta);
    minmaxLB   = betaLB(indB1, indB2) / theta;
    
    fprintf('#%d, min(max(betaLB/theta)): %f, min(max(betaUB/theta)): %f, good: %d/%d\n',...
      zz, minmaxLB, minmaxUB, sum(~tryind), sampleSize);
  end
  
  plot_now;  
 
  if minmaxUB < 1
    k = numel(beta) + 1;
    fprintf('Adding new mubar(%d).\n',k);
    mubar(k,:)  = scmSample(indA1,:);
    muhat{k}    = mubar(k,:);
    LPA{k}      = ThetaVal(indA1,:);
    beta2{k}(1) = 1;
    [beta(k), zvec{k}] = get_infsup(mubar(k,:));
  elseif minmaxLB > opt.epstol,
    break;
  else
    k=indB2;
    indL = numel(beta2{k})+1;
    fprintf('Adding new muhat(%d,%d).\n',k,indL);
    muhat{k}(indL,:)   = scmSample(indB1,:);
    zvec{k}(1:Q, indL) = get_infsup2(muhat{k}(indL,:), mubar(k,:));
    beta2{k}(indL)     = dot(ThetaVal(indB1, :), zvec{k}(:, indL));
    LPA{k}             = [LPA{k}; ThetaVal(indB1,:)];
  end
  save(RB.scm.file,'zvec','beta','gamma','beta2','LPA','mubar','muhat','stats');
end

%% Final bound results
betaSCM = max(bsxfun(@times, betaLB, beta), [], 2);

  function updateUBLB(zvec, beta, gamma, beta2, LPA, whk)
    for n=whk
      time = 0; counter = 0;
      betaUB(:,n) = min(ThetaVal*zvec{n},[],2); % upper bound
      % message
      ntosolve = sum((betaUB(:,n) >= theta) & tryind);
      fprintf('Solving %d linear problems.\n', ntosolve);
      % lower bound
      betaLB(tryind,n) = -Inf; % initial value
      linprogLB = -gamma/beta(n);  linprogA = -LPA{n};
      linprogUB =  gamma/beta(n);  linprogB = -beta2{n};
      parfor ll = 1:sampleSize % lower bound
        % compute the lower bound only if the upper bound is sufficiently
        % large (when traning) or positive (when evaluating) and only for
        % points where sufficient lower bound is not found yet.
        if (betaUB(ll,n) >= train*theta) && tryind(ll)
          tic;
          [~, betaLB(ll,n)] = linprog(ThetaVal(ll,:)', ...
            linprogA,linprogB,[],[],linprogLB,linprogUB,[],...
            struct('Display','none'));
          time = time + toc;
          counter = counter + 1;
        end
      end
      stats.linprog  = stats.linprog  + time; 
      stats.nlinprog = stats.nlinprog + counter; 
    end
  end

  function gamma = get_gamma
    [gamma, time] = deal(zeros(Q,1));    
    parfor qpar=1:Q
      tic;
      Aqfun1 = @(x)(ChT\(S'*(AX{qpar}*(S*(Ch\x)))));
      Aqfun = @(x)(Aqfun1(Aqfun1(x)));
      eigval = bleigifp(Aqfun,1, eigoptLA);
      gamma(qpar) = sqrt(eigval(1));
      time(qpar) = toc;
    end
    stats.eig = [stats.eig; time];
  end

  function zvec = get_infsup2(mu,mubar)
    tic;
    Amu    = assemble_affine(mu,    Theta, AX);
    Amubar = assemble_affine(mubar, Theta, AX);
    Afun1  = @(x)(ChT\(S'*(Amubar*(S*(Ch\x)))));
    Afun2  = @(x)(ChT\(S'*(Amu*(S*(Ch\x)))));
    Afun   = @(x)(1/2*(Afun1(Afun2(x)) + Afun2(Afun1(x))));
    Bfun   = @(x)(Afun1(Afun1(x)));
    [~, eigvec_comp] = bleigifp(Afun, Bfun, 1, eigoptSA);
    [zvec, time] = get_zvec(eigvec_comp(:,1), Afun1(eigvec_comp(:,1)));
    stats.eig = [stats.eig; time+toc];
  end

  function [beta, zvec] = get_infsup(mubar)
    tic;
    Amubar = assemble_affine(mubar, Theta, AX);
    Afun1  = @(x)(ChT\(S'*(Amubar*(S*(Ch\x)))));
    Afun   = @(x)(Afun1(Afun1(x)));
    [eigarray, eigvec_comp] = bleigifp(Afun,1, eigoptSA);
    beta = sqrt(eigarray(1));
    [zvec, time] = get_zvec(eigvec_comp(:,1), Afun1(eigvec_comp(:,1)));
    stats.eig = [stats.eig; time+toc];
  end

  function [zvec, time] = get_zvec(ev1, ev2)
    [zvec, time] = deal(zeros(Q,1));
    ev2norm = dot(ev2,ev2);
    parfor qq=1:Q
      tic;
      ev = (ChT\(S'*(AX{qq}*(S*(Ch\ev1)))));
      zvec(qq) = dot(ev2, ev) / ev2norm;
      time(qq) = toc;
    end
    time = sum(time);
  end

  function plot_now        % plot the current situation
    P = size(mubar, 2);   NS = size(mubar, 1);   colors = jet(NS);
    if P ~= 2, return; end
    goodpts = any(betaLB / theta > opt.epstol, 2);
    [~, father] = max(betaLB / theta,[],2);
    clf;
    scatter(scmSample(goodpts,1), scmSample(goodpts,2), 20, colors(father(goodpts),:));
    hold on;
    for i=1:NS
      for j=1:size(muhat{i},1)
        line([mubar(i,1), muhat{i}(j,1)], [mubar(i,2),muhat{i}(j,2)],'color','black');
      end
    end
    scatter(scmSample(~goodpts,1), scmSample(~goodpts,2),20,[0.7,0.7,0.7],'x');
    scatter(mubar(:,1), mubar(:,2), 100, colors, 'fill');
    for i=1:NS
      scatter(muhat{i}(:,1), muhat{i}(:,2), 50, colors(i,:),'fill');
    end
    drawnow;
  end
end