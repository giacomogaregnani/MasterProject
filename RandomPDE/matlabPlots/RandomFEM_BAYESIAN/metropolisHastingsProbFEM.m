function sampleTot = metropolisHastingsProbFEM(proposal, prior, x, field, f, xObs, observations, sigma, initGuess, N, M)

sample = cell(M, 1);
p = 1;

parfor i = 1 : M
    fprintf('chain number %3d\n', i); 
    X = x;
    for kk = 2 : length(x)-1
        hBar = min(x(kk)-x(kk-1), x(kk+1)-x(kk));
        X(kk) = x(kk) + hBar^p * (rand - 0.5);
    end
    posterior = @(param) prior(param) + likelihoodFEM(param, field, f, X, xObs, observations, sigma);
    sample{i} = metropolisHastings(posterior, proposal, initGuess, N, M);
end

sampleTot = [];
for i = 1 : M
   sampleTot = [sampleTot sample{i}]; 
end

return