%% USER INPUT

% required - normally changed
expName  = 'X';
umicElem = 'p2';
pmicElem = 'p1';
macElem  = 'p2';

% required - not changed often
saveDir = '~/repos/experiments/comp/';
saveMacFile = @(i)(['refsol' num2str(i) '.mat']);

%% DO NOT MODIFY BELOW THIS LINE

% assemble directory names
expDir = [saveDir expName '/'];
expFullName = [umicElem 'vs' umicElem pmicElem];
saveMacDir = [expDir expFullName '/'];

% create non-existent directories
if ~exist(expDir,'dir'),     mkdir(expDir); end
if ~exist(saveMacDir,'dir'), mkdir(saveMacDir); end

% problem definition
macRhs = str2func(['rhs_' expName]);
macVp = struct('elemtype', macElem, 'bc', 0, 'f', macRhs);
S = load(['macMesh_' expName], 'mesh');
macMesh = S.mesh;
macMesh.bdflag = 'neumann';
micMesh = str2func(['micMeshG_' expName]);
micVp = struct('maxDof', 10^6, ...
  'elemtype', umicElem, 'pelemtype', pmicElem, 'bc', [0,0], 'a', 1);


options = struct('saveMacDir', saveMacDir, 'saveMacFile', saveMacFile, 'maxIter', 50, 'maxMicDof', 10^4);

[macSol, macFemspace, macMesh, macVp] = ...
  dsfehmm_uni(macMesh, macVp, micMesh, micVp, options);

%% setting constants (previously computed)
% if strcmp(expName,'c')
%   if strcmp(macElem,'p1') && strcmp(umicElem,'p2') && strcmp(pmicElem,'p1')
%     Cmac = 48.342656;
%     Cmic = 2*0.034661;
%   end
% elseif strcmp(expName,'d')
%   if strcmp(macElem,'p1') && strcmp(umicElem,'p2') && strcmp(pmicElem,'p1')
%     Cmac = 19.018409104833815; 
%     Cmic = 0.081767573391021; 
%   end
% elseif strcmp(expName,'A')
%   if strcmp(macElem,'p1') && strcmp(umicElem,'p2') && strcmp(pmicElem,'p1')
%     Cmac = 29.129247;
%     Cmic = 0.057771; 
%   end
% end
%% RUN
%[Cmic, Cmac, stats] = get_apos_constants(macMesh, macVp, micMesh, micVp, options);

