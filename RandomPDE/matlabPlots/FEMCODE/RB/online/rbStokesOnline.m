function [tensor, errF, stats] = rbStokesOnline(RB, sample)

if ~isfield(RB.online,'estimate'), RB.online.estimate = false; end

%% initialization
sampleSize = size(sample,1);
stats = struct('nrb', sampleSize);

%% LOAD AFFINE DECOMPOSITION (coefficients)
[Theta, ThetaF] = deal([]);  
load(RB.affine.file, 'Theta', 'ThetaF');
Q   = numel(Theta);
QF  = numel(ThetaF);

%% pre-evaluation of Theta fields
ThetaField  = evalTheta(sample, Theta);
ThetaFieldF = evalTheta(sample, ThetaF);

%% SCM FOR INF-SUP LOWER BOUNDS
if RB.online.estimate
  [betaSCM, stats.scm] = natural_SCM(RB, sample, false);
  if any(betaSCM<=0)
    error('negative inf-sup constant lower bound');
  end
end

%% LOAD ONLINE VARIABLES
[OM, ON, OP, OS, OT, N] = deal([]);
load(RB.online.file, 'OM','ON','OP','OS','OT','N');
dim = numel(OM);

%% REDUCED BASIS SOLUTIONS
tic;
redSol = cell(dim,1);
resNorm2 = zeros(sampleSize,dim);
tensor = zeros(sampleSize,dim,dim);
for tr = 1:sampleSize
  T  = ThetaField(tr,:);
  TF = ThetaFieldF(tr,:);
  for i=1:dim
    % assemble the RB system
    Abar = reshape(reshape(reshape(OM{i},[],Q)  * T', [],Q) * T', [N(i),N(i)]);
    Fbar = reshape(reshape(reshape(ON{i},[],QF) * TF',[],Q) * T', [N(i),1]);
    redSol{i} = Abar \ Fbar;
    
    if RB.online.estimate
      resNorm2(tr,i) = ...
        TF * OP{i} * TF' - ...
        2 * (redSol{i}' * Fbar) + ...
        redSol{i}' * Abar * redSol{i};
    end
  end
  for i=1:dim
    for j=1:dim
      tensor(tr,i,j) = ...
        (redSol{j}' * OS{j,i}) * TF' + ...
        (redSol{i}' * OS{i,j}) * TF' - ...
        redSol{i}' * reshape(reshape(OT{i,j},[],Q) * T', [N(i), N(j)]) * redSol{j};        
    end
  end
end
if RB.online.estimate
  resNorm2(resNorm2<0) = 0; % can be negative (up to a round off err)
  errF = sum(resNorm2,2) ./ betaSCM;
else
  errF = [];
end
stats.time = toc;
end