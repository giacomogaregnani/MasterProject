function [macSol, macFemspace, macMesh, macVp] = ...
    dsfehmm_uni(macMesh, macVp, micMeshH, micVp, opt) 
% FEHMM method that saves micro problems

if nargin < 5, opt = struct; end 
if ~isfield(opt, 'maxIter'),   opt.maxIter = 200; end
if ~isfield(opt, 'maxMicDof'), opt.maxMicDof = 10^4; end
if ~isfield(opt, 'maxMacDof'), opt.maxMacDof = 10^6; end
if ~isfield(opt, 'macMark'),   opt.macMark = struct('method', 'L2', 'theta', 0.25); end
if ~isfield(opt, 'micMark'),   opt.micMark = struct('method', 'L2', 'theta', 0.5); end
if ~isfield(opt, 'minMacDiam'), opt.maxMacDof = 10^6; end

global global_options;
GL = global_options; % to be used in parallel workers

%% INITIALIZE
dim = size(macMesh.elem,2) - 1;

macFemspace = get_femspace(macMesh, macVp.elemtype);
macVp.aquad = max(2*(macFemspace.deg-1),1);
[macVp.alambda,macVp.aweight] = quadpts(dim,macVp.aquad);
macVp.fully_discrete = true;
macVp.bc = 0;

NQ  = numel(macVp.aweight);

iMstart = 1;
for i=opt.maxIter:-1:1
  filename = [opt.saveMacDir feval(opt.saveMacFile, i)];
  if exist(filename,'file')
    iMstart = i;
    break;
  end
end

if iMstart > 1
  load(filename,'-regexp', '^(?!opt$).');
  mmp = false(NT, 1, dim, NQ);
else
  NT  = size(macMesh.elem,1);
  mmp = true(NT, 1, dim, NQ); % marked micro problems
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Macro adaptive process
for iM=iMstart : opt.maxIter
  for qp=1:NQ  % for every quadrature point
    for m=1:dim % for every direction
      ind = find(mmp(:,1,m,qp));
      NI = numel(ind);
      xloc = get_rc(macMesh, [], ind, [], macVp.alambda(qp,:));
      average = zeros(NI, dim);
      ndof = zeros(NI, 1);
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %% PARFOR - micro problem solutions (or updates)
      parfor k=1:NI
        setglobal('global_options',GL);
        fprintf('I: %.2d, Dir: %.1d / %.1d, QP: %.1d / %.1d, Pr: %.5d / %.5d\n', ...
          iM, m, dim, qp, NQ, k, NI);
        
        iS = 1;
        micMeshX = feval(micMeshH,xloc(k,:));
        micMeshX.bdflag = 'dirichlet';
        
        micVpX   = micVp;
        micVpX.f = str2func(['fstokesmicro' num2str(m)]);
        [averageX, ndofX] = deal([]);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% STOKES ADAPTIVE SOLUTION ITERATION
        while true
          %% SOLVE
          [micUsol,micUfemspace,micPsol,micPfemspace, ...
            micMeshX, micVpX] = stokes(micMeshX, micVpX);
          
          %% HOMOGENIZED TENSOR
          averageX(iS,:) = get_average(micMeshX, micUfemspace, micUsol);
          ndofX(iS) = dim * micUfemspace.ndof + micPfemspace.ndof;
          
          %% ESTIMATE
          micResX = stokes_residual(micMeshX, micVpX, ...
            micUfemspace, micUsol, micPfemspace, micPsol);
          
          %% IF NUMBER OF DOF REACHED, WRITE TO FILE
          if ndofX(iS) > opt.maxMicDof,
            ndof(k) = ndofX(iS);
            average(k,:) = averageX(iS,:);
            break;
          end
          
          %% MARK AND REFINE MICRO PROBLEM
          me = markelem(micResX, opt.micMark); % MARK
          micMeshX = bisect(micMeshX, me); % REFINE
          iS = iS + 1;
        end %% STOKES ADAPTIVE
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      end %% PARFOR
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      macVp.a(ind,m,:,qp) = average;
      macVp.micNdof(ind,1,m,qp) = ndof;
    end
  end
  %% SYMMETRIZE tensor
  macVp.a = (macVp.a + permute(macVp.a, [1,3,2,4]))/2;
  
  %% SOLVE macro
  [macSol, macFemspace, macMesh, macVp] = poisson(macMesh, macVp);
  
  %% ESTIMATE ERRORS
  etamac = poisson_residual(macMesh, macVp, macFemspace, macSol);
  
  %% MARK AND REFINE
  me = markelem(etamac, opt.macMark);
  diam = simplex_diameter(macMesh);
  marked = me;
  marked(diam < opt.minMacDiam) = 0;
  [newMacMesh, father] = bisect(macMesh, me);
  
  %% SAVE macro iteration
  filename = [opt.saveMacDir feval(opt.saveMacFile, iM)];
  save(filename);
    
  %% CRITERIA TO STOP
  if macFemspace.ndof > opt.maxMacDof || iM >= opt.maxIter,
    break;
  end
  
  %% UPDATE FIELDS IN VP
  NT = size(newMacMesh.elem,1);
  left = false(NT,1);
  left(father(histc(father,unique(father))==1)) = true;
  
  newMacVp = macVp;
  newMacVp.a         = zeros(NT,dim,dim,NQ);
  newMacVp.micNdof   = zeros(NT,1,dim,NQ);
  newMacVp.a(left,:,:,:)          = macVp.a(left,:,:,:);
  newMacVp.micNdof(left,:,:,:)    = macVp.micNdof(left,:,:,:);
  
  macMesh = newMacMesh;
  macVp = newMacVp;
  
  mmp = true(NT, 1, dim, NQ);
  mmp(left,:,:,:) = false;
  clear newMacMesh newMacVp
end
end