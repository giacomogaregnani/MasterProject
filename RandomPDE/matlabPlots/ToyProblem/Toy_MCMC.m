function thetaAll = Toy_MCMC(h, xGrid, prior, likelihood, xObs, nMCMC, nMC, method, uEx, RAM, gamma)

% Initial guess
theta = 1;
oldPrior = prior(theta);

% Initialize structures for MC estimation
if strcmp(method, 'PROB')
    oldLikMC = zeros(1, nMC);
    newLikMC = zeros(1, nMC);
end

% Forward problem
if strcmp(method, 'FEM')
    u = interp1(xGrid, uEx(xGrid, theta), xObs);
    oldLik = likelihood(u);
elseif strcmp(method, 'EX')
    u = uEx(xObs, theta);
    oldLik = likelihood(u);
elseif strcmp(method, 'PROB')
    for j = 1 : nMC
        xGridProb = xGrid + h^2 * rand(size(xGrid));
        u = interp1(xGridProb, uEx(xGridProb, theta), xObs);
        oldLikMC(j) = likelihood(u);
    end
    oldLik = mean(oldLikMC);
end

% Initialize RAM
S = RAMInit(gamma, theta);
accRatio = 0;

% MCMC loop
thetaAll = zeros(1, nMCMC + 1);
thetaAll(1) = theta;

for i = 1 : nMCMC
    
    if mod(i, 5000) == 0
        display(['iteration ' num2str(i), ', current acc. ratio ' num2str(accRatio / i,'%.4f') ...
            ', sigma = ', num2str(max(diag(S)),'%.4f'), ', posterior = ', num2str(oldLik * oldPrior, '%.2f'), ...
            ' theta = ', num2str(thetaAll(i))]);
    end
    
    % Generate new guess (Gaussian proposal)
    w = randn(1);
    theta = thetaAll(i) + S * w;
    
    % Evaluate prior on new guess
    newPrior = prior(theta);
    
    % Solve FWD problem for new guess
    if strcmp(method, 'FEM')
        u = interp1(xGrid, uEx(xGrid, theta), xObs);
        newLik = likelihood(u);

    elseif strcmp(method, 'EX')
        u = uEx(xObs, theta);
        newLik = likelihood(u);

    elseif strcmp(method, 'PROB')
        for j = 1 : nMC
            xGridProb = [0, xGrid(2:end-1) + h^2 * (-0.5 + rand(size(xGrid(2:end-1)))), 1];
            u = interp1(xGrid, uEx(xGridProb, theta), xObs);
            newLikMC(j) = likelihood(u);
            
        end
        newLik = mean(newLikMC);
    end
    
    % Compute alpha
    alpha = newPrior * newLik / (oldPrior * oldLik);
    alpha = min(1, alpha);
    
    % Update the chain
    if alpha > rand(1)
        thetaAll(i+1) = theta;
        oldPrior = newPrior;
        oldLik = newLik;
        accRatio = accRatio + 1;
        thetaOld = theta;
    else
        thetaAll(i+1) = thetaAll(i);
    end
    
    % Update RAM matrix
    if RAM
        S = RAMUpdate(w, alpha, i, 0.234, S);
    end
end
