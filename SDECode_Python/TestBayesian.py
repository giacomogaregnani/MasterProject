from NISDE import *
from ParEst import *
from Hom import *

# Set up the equation
eps = 0.05
sde = MSOrnUhl(eps)
sdeHomo = OrnUhl()
EM = EulerMaruyama(sde)
alpha = np.array([1.0])
sigma = 0.7
IC = 0.0

# Compute the homogenized coefficients
hom_param = compute_homogeneous(sigma, alpha, 2.0 * np.pi, sde.p)
print(hom_param)

np.savetxt("../SDECode/matlabPlots/ParameterEstimation/resultsBayesianOne/ResultsHom.txt", hom_param, fmt='%f')

# Set parameters
betaFilter = 1
gammaTime = -np.log(2000) / np.log(eps)
beta = 3.0
T = round(eps ** (-1.0 * gammaTime))
h = eps ** beta
N = int(round(T/h))
print('T = ', T, 'h = ', h, 'N = ', N)
delta = eps ** 0.5

# Generate solution
y, dw = EM.solve(IC, N, T, alpha, sigma, seed=2020)

# Filter and estimate parameter
filt, kt = filter_trajectory(y, delta, T, beta=betaFilter, type=1)
par_est_filt = ParEst(filt, sdeHomo.grad_v, h, is_vect=True)
meanPrior = 0.
varPrior = 1.
postMean, postCov = par_est_filt.drift_mixed_bayesian_all(hom_param[-1], meanPrior, varPrior, data=y)
aHat = par_est_filt.drift_mixed(data=y)

np.savetxt("../SDECode/matlabPlots/ParameterEstimation/resultsBayesianOne/ResultsFilterT2_" + str(int(T)) + "MLE.txt", [aHat], fmt='%f')
np.savetxt("../SDECode/matlabPlots/ParameterEstimation/resultsBayesianOne/ResultsFilterT2_" + str(int(T)) + "Mean.txt", [postMean], fmt='%f')
np.savetxt("../SDECode/matlabPlots/ParameterEstimation/resultsBayesianOne/ResultsFilterT2_" + str(int(T)) + "Cova.txt", [postCov], fmt='%f')

