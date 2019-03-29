function [macSol, macFemspace, macMesh, macVp, stats] = ...
    dsfehmm(macMesh, macVp, micMeshH, micVp, opt) 
% FEHMM method that saves micro problems

if nargin < 5, opt = struct; end 
if ~isfield(opt, 'maxIter'),   opt.maxIter = 200; end
if ~isfield(opt, 'maxInIter'), opt.maxInIter = 200; end
if ~isfield(opt, 'maxMicDof'), opt.maxMicDof = 10^6; end
if ~isfield(opt, 'minMicDof'), opt.minMicDof = 0; end
if ~isfield(opt, 'maxMacDof'), opt.maxMacDof = 10^6; end
if ~isfield(opt, 'muMac'),     opt.muMac = 1; end
if ~isfield(opt, 'muMic'),     opt.muMic = 1; end
if ~isfield(opt, 'muQ'),       opt.muQ = 1; end
if ~isfield(opt, 'saveMac'),   opt.saveMac = true; end
if ~isfield(opt, 'saveMic'),   opt.saveMic = true; end
if ~isfield(opt, 'macMark'),   opt.macMark = struct('method', 'L2', 'theta', 0.25); end
if ~isfield(opt, 'micMark'),   opt.micMark = struct('method', 'L2', 'theta', 0.5); end
if opt.saveMac && ~isfield(opt, 'saveMacFile'),
  error('You must specify where macro solutions are to be stored.');
end
if opt.saveMic && ~isfield(opt, 'saveMicFile')
  error('You must specify where macro solutions are to be stored.');
end

global global_options;
GL = global_options; % to be used in parallel workers

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
fmingradp = ones(NT,1);     % | f - grad p |_L2(K)
etamac = Inf(NT,1);         % macro residual
mmp = true(NT, 1, dim, NQ); % marked micro problems
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Macro adaptive process
for iM=1:opt.maxIter
  stats{iM}.micro = 0;
  stats{iM}.macro = 0;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% Inner iterations - to make error(mic) < error(mac)
  for iI = 1:opt.maxInIter
    for qp=1:NQ  % for every quadrature point
      for m=1:dim % for every direction
        ind = find(mmp(:,1,m,qp));
        NI = numel(ind);
        xloc = get_rc(macMesh, [], ind, [], macVp.alambda(qp,:));
        micStop = 1/dim * opt.muMac^2 / (opt.muQ * opt.muMic)^2 * etamac(ind) ./ fmingradp(ind);
        
        etaStokes = zeros(NI, 1);
        average = zeros(NI, dim);
        ndof = zeros(NI, 1);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% PARFOR - micro problem solutions (or updates)
        time = 0;
        parfor k=1:NI
          setglobal('global_options',GL);
          fprintf('I: %.2d, In. I: %.2d, Dir: %.1d / %.1d, QP: %.1d / %.1d, Pr: %.5d / %.5d\n', ...
            iM, iI, m, dim, qp, NQ, k, NI);
          saveMicFile = [opt.saveMicDir, feval(opt.saveMicFile, xloc(k,:), m, qp)];
          
          try % load
            [etaStokesX, averageX, ndofX] = load_from_file(saveMicFile, 'etaStokes', 'average', 'ndof');
            iS = find(((etaStokesX < micStop(k)) | (ndofX >= opt.maxMicDof)) &...
              ndofX >= opt.minMicDof, 1);
            if isempty(iS) && (ndofX(end) >= opt.maxMicDof)
              iS = numel(ndofX);
            end
            if ~isempty(iS)
              ndof(k) = ndofX(iS);
              etaStokes(k) = etaStokesX(iS);
              average(k,:) = averageX(iS,:);
              continue;
            else
              [micMeshX, micResX, micVpX] = load_from_file(saveMicFile, 'micMesh', 'micRes', 'micVp');
              me = markelem(micResX); % MARK
              micMeshX = bisect(micMeshX, me); % REFINE
              iS = numel(ndofX) + 1;
            end
          catch % FIRST TIME SOLVING THIS MICRO PROBLEM
            iS = 1;
            micMeshX = feval(micMeshH,xloc(k,:));
            micMeshX.bdflag = 'dirichlet';
            
            micVpX   = micVp;
            micVpX.f = str2func(['fstokesmicro' num2str(m)]);
            [etaStokesX, averageX, ndofX] = deal([]);
          end
          
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          %% STOKES ADAPTIVE SOLUTION ITERATION
          tic;
          while true
            %% SOLVE
            [micUsol,micUfemspace,micPsol,micPfemspace, ...
              micMeshX, micVpX] = stokes(micMeshX, micVpX);
            
            %% MINI ELEMENTS - remove bubbles
            if strcmp(micUfemspace.elemtype,'p1b')
              Nmic = size(micMeshX.node,1);
              micUsol = micUsol(Nmic,:,:);
              micUfemspace = get_femspace(micMeshX,'p1');
            end
            
            %% HOMOGENIZED TENSOR
            averageX(iS,:) = get_average(micMeshX, micUfemspace, micUsol);
            ndofX(iS) = dim * micUfemspace.ndof + micPfemspace.ndof;
            
            %% ESTIMATE
            micResX = stokes_residual(micMeshX, micVpX, ...
              micUfemspace, micUsol, micPfemspace, micPsol);
            etaStokesX(iS) = sum(micResX, 1);
            
            %% IF PRECISION OR DOF REACHED, WRITE TO FILE
            if (etaStokesX(iS) < micStop(k) && ndofX(iS) >= opt.minMicDof) || ...
                ndofX(iS) > opt.maxMicDof,
              time = time + toc;
              if opt.saveMic
                save_to_file(saveMicFile, ...
                  'micMesh', micMeshX,'micRes', micResX,...
                  'micVp', micVpX, 'etaStokes', etaStokesX,...
                  'average',averageX,'ndof', ndofX);
                ndof(k) = ndofX(iS);
              end
              etaStokes(k) = etaStokesX(iS);
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
        stats{iM}.micro = stats{iM}.micro + time; 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        macVp.a(ind,m,:,qp) = average;
        macVp.etaStokes(ind,1,m,qp) = etaStokes;
        macVp.micNdof(ind,1,m,qp) = ndof;
      end
    end
    %% SYMMETRIZE tensor
    macVp.a = (macVp.a + permute(macVp.a, [1,3,2,4]))/2;    
    
    %% SOLVE macro
    tic;
    [macSol, macFemspace, macMesh, macVp] = poisson(macMesh, macVp);
    
    %% ESTIMATE ERRORS
    etamac = poisson_residual(macMesh, macVp, macFemspace, macSol);
    fmingradp = get_fmingradp(macMesh, macVp, macFemspace, macSol);
    etamic = fmingradp .* max(sum(macVp.etaStokes,3),[],4);
    
    %% Mark micro problems
    % where the micro residual is too big and the maximal number of degrees
    % of freedom was not yet reached
    mmp = bsxfun(@gt, macVp.etaStokes, ...
      1/dim * opt.muMac^2 / (opt.muQ * opt.muMic)^2 ...
      * etamac ./ fmingradp) & (macVp.micNdof < opt.maxMicDof);
    stats{iM}.macro = stats{iM}.macro + toc;
    if ~any(mmp(:)), break; end
  end
  
%% MARK AND REFINE
  tic; 
  me = markelem(etamac, opt.macMark);
  [newMacMesh, father] = bisect(macMesh, me);
  stats{iM}.macro = stats{iM}.macro + toc;
  
  %% SAVE macro iteration
  if opt.saveMac
    filename = [opt.saveMacDir feval(opt.saveMacFile, iM)];
    save(filename);
  end
  
  %% CRITERIA TO STOP
  if macFemspace.ndof > opt.maxMacDof || iM >= opt.maxIter, 
    break; 
  end
  
  %% UPDATE FIELDS IN VP
  tic;
  NT = size(newMacMesh.elem,1);
  left = false(NT,1);
  left(father(histc(father,unique(father))==1)) = true;
  
  newMacVp = macVp;
  newMacVp.a         = zeros(NT,dim,dim,NQ);
  newMacVp.etaStokes = zeros(NT,1,dim,NQ);
  newMacVp.micNdof   = zeros(NT,1,dim,NQ);
  newMacVp.a(left,:,:,:)          = macVp.a(left,:,:,:);
  newMacVp.etaStokes(left,:,:,:)  = macVp.etaStokes(left,:,:,:);
  newMacVp.micNdof(left,:,:,:)    = macVp.micNdof(left,:,:,:);
  
  etamic = etamic(father);
  etamac = etamac(father); % preliminary values for the next iteration
  fmingradp = fmingradp(father);
  
  macMesh = newMacMesh;
  macVp = newMacVp;
  
  mmp = true(NT, 1, dim, NQ);
  mmp(left,:,:,:) = false;  
  clear newMacMesh newMacVp
  stats{iM}.macro = stats{iM}.macro + toc;
end
end