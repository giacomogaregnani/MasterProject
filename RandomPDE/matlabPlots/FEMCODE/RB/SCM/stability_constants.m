function betaSCM = stability_constants(RB, sample)

opt = RB.scm;
if ~isfield(opt,'eigtol'), opt.eigtol = 1e-4; end
sampleSize = size(sample,1);

%% LOAD + INIT + SETTINGS
[Ch, S, AX, Theta] = deal([]);
load(RB.affine.file, 'Ch','S','AX','Theta');
ChT = Ch';
NX = size(AX{1},1);
eigoptSA = struct('tol', opt.eigtol, 'size',NX,'sigma','SA');

%% SCM MAIN ITERATION
betaSCM = zeros(sampleSize,1);
for k=1:sampleSize % this will end much sooner
    betaSCM(k) = get_infsup(sample(k,:));
end

  function beta = get_infsup(mubar)
    Amubar = assemble_affine(mubar, Theta, AX);
    Afun1  = @(x)(ChT\(S'*(Amubar*(S*(Ch\x)))));
    Afun   = @(x)(Afun1(Afun1(x)));
    eigarray = bleigifp(Afun,1, eigoptSA);
    beta = sqrt(eigarray(1));    
  end
end