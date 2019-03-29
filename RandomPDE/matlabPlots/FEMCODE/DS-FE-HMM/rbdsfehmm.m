function [macSol, macFemspace, macMesh, macVp, stats] = ...
    rbdsfehmm(macMesh, macVp, RB, opt) 
  
if nargin < 4, opt = struct; end
if ~isfield(opt, 'maxIter'),   opt.maxIter = 200; end
if ~isfield(opt, 'maxMacDof'), opt.maxMacDof = 10^6; end
if ~isfield(opt, 'macMark'),   opt.macMark = struct('method', 'L2', 'theta', 0.25); end
if ~isfield(opt, 'saveMac'),   opt.saveMac = true; end
  
if ~isfield(RB.online,'estimate'), RB.online.estimate = false; end
  
  %% INITIALIZE
  dim = size(macMesh.elem,2) - 1;
  stats = cell(1,1);

macFemspace = get_femspace(macMesh, macVp.elemtype);
macVp.aquad = max(2*(macFemspace.deg-1),1);
[macVp.alambda,macVp.aweight] = quadpts(dim,macVp.aquad);
macVp.fully_discrete = true;
macVp.bc = 0;

NQ  = numel(macVp.aweight);
NT  = size(macMesh.elem,1);
newelem = (1:NT)';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Macro adaptive process
for iM=1:opt.maxIter
  NI = numel(newelem);
  xloc = zeros(NQ*NI,dim);
  for qp=1:NQ  % for every quadrature point
      xloc(1+(qp-1)*NI:qp*NI,:) = ...
        get_rc(macMesh, [], newelem, [], macVp.alambda(qp,:));
  end
  fprintf('I: %.2d, mic problems: %d, ', iM, size(xloc,1));
  [tensor, errF, stats{iM}.online] = rbStokesOnline(RB, macVp.x2par(xloc));
  for qp=1:NQ
    macVp.a(newelem,:,:,qp) = tensor(1+(qp-1)*NI:qp*NI,:,:);
    if RB.online.estimate
      macVp.errF(newelem,qp) = errF(1+(qp-1)*NI:qp*NI);
    end
  end  
  
  %% SYMMETRIZE tensor
  macVp.a = (macVp.a + permute(macVp.a, [1,3,2,4]))/2;
  
  %% SOLVE macro
  tic;
  [macSol, macFemspace, macMesh, macVp] = poisson(macMesh, macVp);
  fprintf('dof: %d/%d\n', macFemspace.ndof, opt.maxMacDof);
  
  %% ESTIMATE ERRORS
  etamac = poisson_residual(macMesh, macVp, macFemspace, macSol);
  if RB.online.estimate
    fmingradp = get_fmingradp(macMesh, macVp, macFemspace, macSol);
    etamic = fmingradp .* max(macVp.errF, [], 2);
    if any(etamic > etamac)
      fprintf('There are %d elements with etamic > etamac\n', ...
        sum(etamic>etamac));
    end
  end
   
  %% MARK AND REFINE
  me = markelem(etamac, opt.macMark);
  [newMacMesh, father] = bisect(macMesh, me);
  stats{iM}.solve = toc;
  
  %% SAVE macro iteration
  if opt.saveMac
    filename = [opt.saveMacDir feval(opt.saveMacFile, iM)];
    save(filename);
  end
  
  %% CRITERIA TO STOP
  if macFemspace.ndof > opt.maxMacDof || iM >= opt.maxIter, break; end
  
  %% UPDATE FIELDS IN VP
  NT = size(newMacMesh.elem,1);
  left = false(NT,1);
  left(father(histc(father,unique(father))==1)) = true;
  
  newMacVp = macVp;
  newMacVp.a    = zeros(NT,dim,dim,NQ);
  newMacVp.a(left,:,:,:) = macVp.a(left,:,:,:);
  if RB.online.estimate
    newMacVp.errF = zeros(NT,NQ);
    newMacVp.errF(left,:)  = macVp.errF(left,:);
  end
  
  macMesh = newMacMesh;
  macVp = newMacVp;  
  
  newelem = find(~left);
  clear newMacMesh newMacVp
end
end