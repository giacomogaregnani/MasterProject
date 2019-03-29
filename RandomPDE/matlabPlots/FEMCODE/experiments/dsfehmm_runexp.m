%% USER INPUT

% required - normally changed
expName  = 'X';
umicElem = 'p2';
pmicElem = 'p1';
macElem  = 'p1';

% required - not changed often
saveDir = '~/repos/experiments/comp/';
CQ = 2;
saveMicFile = @generate_filename;
saveMacFile = @(i)(['iter' num2str(i) '.mat']);

%% DO NOT MODIFY BELOW THIS LINE

% assemble directory names
expDir = [saveDir expName '/'];
expFullName = [macElem 'vs' umicElem pmicElem];
saveMacDir = [expDir expFullName '/'];
saveMicDir = [expDir umicElem pmicElem '/'];
constantsFile = [saveMacDir 'constants.mat'];

% create non-existent directories
if ~exist(expDir,'dir'),     mkdir(expDir); end
if ~exist(saveMacDir,'dir'), mkdir(saveMacDir); end
if ~exist(saveMicDir,'dir'), mkdir(saveMicDir); end

% problem definition
macRhs = str2func(['rhs_' expName]);
macVp = struct('elemtype', macElem, 'bc', 0, 'f', macRhs);
S = load(['macMesh_' expName], 'mesh');
macMesh = S.mesh;
macMesh.bdflag = 'neumann';
micMesh = str2func(['micMesh_' expName]);
micVp = struct('maxDof', 10^6, ...
  'elemtype', umicElem, 'pelemtype', pmicElem, 'bc', [0,0], 'a', 1);


try 
  load(constantsFile,'Cmac','Cmic')
catch
  options = struct;
  [Cmic, Cmac, stats] = get_apos_constants(macMesh, macVp, micMesh, micVp, options);
  save(constantsFile,'Cmic','Cmac','stats');
end

options = struct('saveMicDir', saveMicDir, 'saveMicFile', saveMicFile, ...
  'saveMacDir', saveMacDir, 'saveMacFile', saveMacFile, ...
  'muMac', Cmac, 'muMic', Cmic, 'muQ', CQ, 'maxIter', 50);

[macSol, macFemspace, macMesh, macVp, stats] = ...
  dsfehmm(macMesh, macVp, micMesh, micVp, options);

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

