function [tensor, errF, sol] = rbStokesOnlineFull(RB, sample, M)
% M = number of basis functions to use, vector of size dim.

error('This function is outdated - see changes in rbStokesOnline.m')

%% initialization
sampleSize = size(sample,1);


%% LOAD AFFINE DECOMPOSITION (coefficients)
[Theta, ThetaF] = deal([]);  
load(RB.affine.file, 'Theta', 'ThetaF');
Q   = numel(Theta);
QF  = numel(ThetaF);

%% pre-evaluation of Theta fields
ThetaField  = evalTheta(sample, Theta);
ThetaFieldF = evalTheta(sample, ThetaF);

%% SCM FOR INF-SUP LOWER BOUNDS
betaSCM = natural_SCM(RB, sample, false);
if any(betaSCM<=0)
  error('negative inf-sup constant lower bound');
end


%% LOAD ONLINE VARIABLES
[OM, ON, OP, OQ, OR, OS, OT, OU, N, RED] = deal([]);
load(RB.online.file, 'OM','ON','OP','OQ','OR','OS','OT','OU','N');
load(RB.offline.file, 'RED');
dim = numel(OM);

if nargin < 3
  M = N;
end

%% REDUCED BASIS SOLUTIONS
tic;
redSol = cell(dim,1);
resNorm2 = zeros(sampleSize,dim);
tensor = zeros(sampleSize,dim,dim);
sol = cell(sampleSize,1);
for tr = 1:sampleSize
  T  = ThetaField(tr,:);
  TF = ThetaFieldF(tr,:);
  for i=1:dim
    % assemble the RB system
    Abar = reshape(reshape(reshape(OM{i}(1:M(i),1:M(i),:,:),[],Q)  * T', [],Q) * T', [M(i),M(i)]);
    Fbar = reshape(reshape(reshape(ON{i}(1:M(i),:,:),[],QF) * TF',[],Q) * T', [M(i),1]);
    redSol{i} = Abar \ Fbar;
    
    resNorm2(tr,i) = ...
      TF * OP{i} * TF' + ...
      2*((redSol{i}' * reshape(reshape(OQ{i}(1:M(i),:,:),[],QF) * TF',M(i),Q)) * T') + ...
      redSol{i}' * reshape(reshape(reshape(OR{i}(1:M(i), 1:M(i),:,:),[],Q) * T', [], Q) * T', M(i), M(i)) * redSol{i};
        
    sol{tr,i} = RED{i}(:,1:M(i)) * redSol{i};    
  end
  
  
  for i=1:dim
    for j=1:dim
      tensor(tr,i,j) = ...
        TF * (OU{i,j}(:,1:M(j)) * redSol{j}) + ...
        ((redSol{i}' * OS{i,j}(1:M(i),:)) * TF') + ...
        redSol{i}' * reshape(reshape(OT{i,j}(1:M(i),1:M(j),:),[],Q) * T', [M(i), M(j)]) * redSol{j};
    end
  end
end

resNorm2(resNorm2<0) = 0; % can be negative (up to a round off err)
errF = sum(resNorm2,2) ./ betaSCM;
end