function [Cmic, Cmac, stats] = ...
  get_apos_constants(macMesh, macVp, micMeshH, micVp, opt) 
%GET_APOS_CONSTANTS computes the a posteriori constants by brute-force

if nargin < 5, opt = struct; end 
if ~isfield(opt, 'macRef1'),   opt.macRef1 = 0; end
if ~isfield(opt, 'macRef2'),   opt.macRef2 = 2; end
if ~isfield(opt, 'macElem'),   opt.macElem = 'p2'; end
if ~isfield(opt, 'micRef1'),   opt.micRef1 = 0; end
if ~isfield(opt, 'micRef2'),   opt.micRef2 = 3; end
if ~isfield(opt, 'micElem'),   opt.micElem = 'p2'; end
if ~isfield(opt, 'micPelem'),  opt.micPelem = 'p1'; end

global global_options;
GL = global_options; % to be used in parallel workers

%% INITIALIZE 
dim = size(macMesh.elem,2) - 1;
[macMeshX, macVpX,macSolX, macFemspaceX, fmingradp, etamac, etamic] = deal(cell(3,1));
stats.macro = 0;
stats.micro = 0;

for pr = 1:3  
  tic;
  macMeshX{pr} = macMesh;
  for k=1 : opt.macRef1 
    macMeshX{pr} = uniformrefine(macMeshX{pr});
  end
  if pr == 3
    NT = size(macMeshX{pr}.elem,1);  
    father = (1:NT)';
    for k = opt.macRef1+1 : opt.macRef2
      [macMeshX{pr}, father] = uniformrefine(macMeshX{pr}, father);
    end
  end

  macVpX{pr} = macVp;
  if pr == 3
    macVpX{pr}.elemtype = opt.macElem;
  end

  macFemspaceX{pr} = get_femspace(macMeshX{pr}, macVpX{pr}.elemtype);
  macVpX{pr}.aquad = max(2*(macFemspaceX{pr}.deg-1),1);
  [macVpX{pr}.alambda, macVpX{pr}.aweight] = quadpts(dim,macVpX{pr}.aquad);
  macVpX{pr}.fully_discrete = true;
  macVpX{pr}.bc = 0;

  NQ  = numel(macVpX{pr}.aweight);
  NT  = size(macMeshX{pr}.elem,1);  
  stats.macro = stats.macro + toc;
  
  for qp=1:NQ  % for every quadrature point
    xloc = get_rc(macMeshX{pr}, [], 'all', [], macVpX{pr}.alambda(qp,:));
    for m=1:dim % for every direction
      average = zeros(NT,dim);
      etaStokes = zeros(NT,1);
      time = 0;
      parfor i=1:NT
        tic;
        setglobal('global_options',GL);
        fprintf('Dir: %.1d / %.1d, QP: %.1d / %.1d, Pr: %.5d / %.5d\n', ...
          m, dim, qp, NQ, i, NT);
        
        %% DEFINE MICRO MESH AND VARIATIONAL PROBLEM
        micMeshX = feval(micMeshH, xloc(i,:));
        for k=1 : (pr<2) * opt.micRef1 + (pr >= 2) * opt.micRef2
          micMeshX = uniformrefine(micMeshX);
        end     
        micMeshX.bdflag = 'dirichlet';
        
        micVpX   = micVp;
        micVpX.f = str2func(['fstokesmicro' num2str(m)]);
        if pr >= 2 
          micVpX.elemtype = opt.micElem;
          micVpX.pelemtype = opt.micPelem;
        end
        
        %% SOLVE
        [micUsol, micUfemspace, micPsol, micPfemspace, ...
          micMeshX, micVpX] = stokes(micMeshX, micVpX);
        
        %% MINI ELEMENTS - remove bubbles
        if strcmp(micUfemspace.elemtype,'p1b')
          Nmic = size(micMeshX.node,1);
          micUsol = micUsol(Nmic,:,:);
          micUfemspace = get_femspace(micMeshX,'p1');
        end
        %% HOMOGENIZED TENSOR
        average(i,:) = get_average(micMeshX, micUfemspace, micUsol);
        
        %% ESTIMATE
        micRes = stokes_residual(micMeshX, micVpX, ...
          micUfemspace, micUsol, micPfemspace, micPsol);
        etaStokes(i) = sum(micRes, 1);
        time = time + toc;
      end
      stats.micro = stats.micro + time;
      macVpX{pr}.a(1:NT,m,1:dim,qp) = reshape(average, [NT, 1, dim, 1]);
      macVpX{pr}.etaStokes(1:NT,1,m,qp) = etaStokes;
    end
  end
  
  %% SYMMETRIZE tensor
  macVpX{pr}.a = (macVpX{pr}.a + permute(macVpX{pr}.a, [1,3,2,4]))/2;
  
  %% SOLVE macro and ESTIMATE ERRORS
  tic;
  [macSolX{pr}, macFemspaceX{pr}, macMeshX{pr}, macVpX{pr}] = poisson(macMeshX{pr}, macVpX{pr});
  etamac{pr} = poisson_residual(macMeshX{pr}, macVpX{pr}, macFemspaceX{pr}, macSolX{pr});
  fmingradp{pr} = get_fmingradp(macMeshX{pr}, macVpX{pr}, macFemspaceX{pr}, macSolX{pr});
  etamic{pr} = fmingradp{pr} .* max(sum(macVpX{pr}.etaStokes,3),[],4);
  stats.macro = stats.macro + toc;
end

tic;
NT = size(macMeshX{1}.elem,1);
[ ~, ~, h1semierrmic] = ...
	get_H1error(macMeshX{1}, macFemspaceX{1}, macSolX{1}, ...
	macMeshX{2}, macFemspaceX{2}, macSolX{2}, (1:NT)');
[ ~, ~, h1semierrmac] = ...
	get_H1error(macMeshX{2}, macFemspaceX{2}, macSolX{2}, ...
	macMeshX{3}, macFemspaceX{3}, macSolX{3}, father);
resmacest = sqrt(sum(etamac{1},1));
resmicest = sqrt(sum(etamic{1},1));
stats.macro = stats.macro + toc;

%% OUTPUTS
Cmac = h1semierrmac/resmacest;
Cmic = h1semierrmic/resmicest;
fprintf(1,'Cmac = %f\n',Cmac);
fprintf(1,'Cmic = %f\n',Cmic);
end

