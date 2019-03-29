function MCMCEst = plotMCMCResultsPDE(samples, uMean, V, sqrtD, MESH)

% Results
MCMCEst = mean(samples, 2);
MCMCDiffusion = exp(KL(uMean, V, sqrtD, MCMCEst));

MCMCDiffusion2 = zeros(length(MCMCDiffusion), size(samples, 2));
for i = 1 : size(samples, 2)
   MCMCDiffusion2(:, i) = KL(uMean, V, sqrtD, samples(:, i));
end
MCMCDiffusion2Mean = mean(exp(MCMCDiffusion2), 2);
MCMCDiffusion2Var = var(exp(MCMCDiffusion2), 0, 2);
clear MCMCDiffusion2

FEMplot2D(MESH, MCMCDiffusion)
title('P(E[\theta])')
FEMplot2D(MESH, MCMCDiffusion2Mean)
title('E[P(\theta)]')
FEMplot2D(MESH, MCMCDiffusion2Var)
title('Var[P(\theta)]')