function sample = run_bayesian(Nsample, noisy_observations, std_noise, locations, log_prior, mesh, prob, errEst, nErrEst, nFine, nIter)

if nargin == 7
    errEst = false;
end

if errEst
    
    sample = [];
    for j = 1 : nIter
        if j == 1
            [mErr, sErr] = estimateError(mesh, locations, nErrEst, nFine);
        else
            [mErr, sErr] = estimateError(mesh, locations, nErrEst, nFine, newSample(:, randsample(1:size(newSample, 2), nErrEst)));
        end
        newNoise = std_noise^2 * eye(size(sErr)) + sErr;
        log_likelihood = @(p) log_likelihood_homog(p, mesh, locations, noisy_observations(:) - mErr, newNoise);
        noisy = false;
        log_posterior = @(p) log_prior(p) + log_likelihood(p);
        newSample = Metropolis_Hastings([log(7); 0.3], log_posterior, 0.002, Nsample, noisy);
        sample = [sample, newSample];
    end
    
elseif prob
    nMC = 20;
    log_likelihood = @(p) log_likelihood_homog_prob(p, mesh, locations, noisy_observations(:), std_noise, nMC);
    noisy = true;
    log_posterior = @(p) log_prior(p) + log_likelihood(p);
    sample = Metropolis_Hastings([log(7); 0.3], log_posterior, 0.002, Nsample, noisy);
    
else
    log_likelihood = @(p) log_likelihood_homog(p, mesh, locations, noisy_observations(:), std_noise);
    noisy = false;
    log_posterior = @(p) log_prior(p) + log_likelihood(p);
    sample = Metropolis_Hastings([log(7); 0.3], log_posterior, 0.002, Nsample, noisy);
    
end

return