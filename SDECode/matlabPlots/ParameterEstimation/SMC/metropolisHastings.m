function sample = metropolisHastings(posterior, proposal, initGuess, N)
% Only for symmetric proposals
% Posterior is a log-posterior

nParam = length(initGuess);
sample = zeros(nParam, N);
sample(:, 1) = initGuess;
posteriorOld = posterior(initGuess);

accRatio = 0;

for i = 1 : N
        
    thetaOld = sample(:, i);
    
    thetaNew = proposal(thetaOld);
    posteriorNew = posterior(thetaNew);
   
    alpha = min(0, posteriorNew - posteriorOld);
    
    if alpha > log(rand)
        sample(:, i+1) = thetaNew;
        accRatio = accRatio + 1;
        posteriorOld = posteriorNew;
    else
        sample(:, i+1) = thetaOld;
    end
    
    if ~mod(i, 50)
        fprintf('iterations completed: %d, acceptance ratio: %.2f, posterior: %.2f \n', i, accRatio/i, posteriorOld)
    end
    
end

cutOff = N / 10;
sample = sample(:, cutOff+1 : end);

end