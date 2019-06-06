%% USER INPUT
expName = 'X';
umicElem = 'p2';
pmicElem = 'p1';
getFineAlg = 'ur2';
x2parAlg = 'x2par';
useCholesky = true;

macElem  = 'p1';

%% RB structure
% RB.expName     code name for the experiment
% RB.expFolder   folder where everything related to the experiment RB.expName is.
% RB.elemtype    elements used for velocity
% RB.pelemtype   elements used for pressure
% RB.getFineAlg  code name for the algorithm that
% RB.affine         structure with info regarding the affine decomposition
% RB.affine.getfine    function handle that takes the coarse mesh and generates the fine mesh
% RB.affine.file       path to mat-file with the affine decomposition
% RB.scm         structure with info regarding the SCM
% RB.scm.file    path to mat-file with scm info
% RB.scm.theta   theta-parameter from the SCM, lies in (0,1)
% RB.scm.epstol  tolerance of the SCM, lies in (0,1)
% RB.scm.eigtol  eigenvalue tolerance for SCM, default 1e-3

global global_options;
compFolder = global_options.comp;

%% creating RB structure
RB = struct('expName', expName, 'elemtype', umicElem, ...
  'pelemtype', pmicElem, 'getFineAlg', getFineAlg,'chol',useCholesky);
RB.problemString = ['rb' RB.elemtype RB.pelemtype RB.getFineAlg];
RB.expFolder = [compFolder RB.expName '/' RB.problemString '/'];
RB.getfine = str2func(['get_fine_' RB.getFineAlg '_' RB.expName]);
RB.reduce  = str2func(['reduce_' RB.expName]);

if ~exist(RB.expFolder,'dir'), mkdir(RB.expFolder); end

saveMacDir = [RB.expFolder x2parAlg '/'];
if ~exist(saveMacDir,'dir'), mkdir(saveMacDir); end
saveMacFile = @(i)(['iter' num2str(i) '.mat']);

%% initialize
RB = feval(['rbinit_' RB.expName], RB);
[RB.affine, RB.scm, RB.offline, RB.online] = deal(struct);
RB.affine.file  = [RB.expFolder 'affine'  RB.expName '.mat'];
RB.scm.file     = [RB.expFolder 'scm'     RB.expName '.mat'];
RB.offline.file = [RB.expFolder 'offline' RB.expName '.mat'];
RB.offline.save = false;
RB.online.file  = [RB.expFolder 'online'  RB.expName '.mat'];

%% SCM options
RB.scm = catstruct(RB.scm, struct('theta', 0.5, 'epstol', 0.5, 'eigtol', 1e-5));

%% RB Greedy settings
RB.offline = catstruct(RB.offline, struct('epstol', 1e-6, 'fakescm', false));

if ~exist(RB.affine.file, 'file'),
  decompose(RB);
end
if ~exist(RB.scm.file, 'file')
  sampleSettings = struct('Npoints', 201^2, 'regular', true);
  scmSample = generate_sample(RB.param, sampleSettings);
  [betaSCM, statsSCMoffline] = natural_SCM(RB, scmSample, true);
end
if ~exist(RB.online.file, 'file')  
  sampleSettings = struct('Npoints', 101^2, 'regular', false);
  rbSample = generate_sample(RB.param, sampleSettings);
  statsRBoffline = rbStokesOffline(RB, rbSample);
end

%% multiscale problem
macRhs = str2func(['rhs_' RB.expName]);
x2par = str2func([x2parAlg '_' RB.expName]);
macVp = struct('elemtype', macElem, 'bc', 0, 'f', macRhs, 'x2par', x2par);
S = load(['macMesh_' RB.expName], 'mesh');
macMesh = S.mesh;
macMesh.bdflag = 'neumann';

RB.online.estimate = false;
opt = struct('maxMacDof',10^4,'saveMac', true, 'saveMacDir', saveMacDir, 'saveMacFile', saveMacFile);
[macSol, macFemspace, macMesh, macVp, stats] = ...
  rbdsfehmm(macMesh, macVp, RB, opt);