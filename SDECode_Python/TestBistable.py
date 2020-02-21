from NISDE import *
from ParEst import *
from Hom import *

# Set up the equation
eps = 0.05
sde = MSBistable2(eps)
sdeHomo = Bistable2()
EM = EulerMaruyama(sde)
alpha = np.array([1.0, 2.0])
sigma = 0.7
IC = 0.0

# Compute the homogenized coefficients
hom_param = compute_homogeneous(sigma, alpha, 2.0 * np.pi, sde.p)
print(hom_param)

np.savetxt("../SDECode/matlabPlots/ParameterEstimation/resultsTwo/ResultsHom.txt", hom_param, fmt='%f')

# Set parameters
betaFilter = 1
gammaTime = -np.log(100) / np.log(eps)
beta = 3.0
T = round(eps ** (-1.0 * gammaTime))
h = eps ** beta
N = int(round(T/h))
print('T = ', T, 'h = ', h, 'N = ', N)
delta = eps ** 1.0

# Generate solution
y, dw = EM.solve(IC, N, T, alpha, sigma, seed=2019)

# Filter and estimate parameter
filt, kt = filter_trajectory(y, delta, T, beta=betaFilter, type=1)
par_est_filt = ParEst(filt, sdeHomo.grad_v, h, is_vect=True)
meanPrior = np.zeros(2)
covPrior = np.eye(2)
postMean, postCov = par_est_filt.drift_multi_mixed_bayesian(2, hom_param[-1], meanPrior, covPrior, data=y, tilde=True)
aHat = par_est_filt.drift_multi_mixed(2, data=y, tilde=True)

np.savetxt("../SDECode/matlabPlots/ParameterEstimation/resultsTwo/ResultsFilterT_" + str(int(T)) + "MLE.txt", aHat, fmt='%f')
np.savetxt("../SDECode/matlabPlots/ParameterEstimation/resultsTwo/ResultsFilterT_" + str(int(T)) + "Mean.txt", postMean, fmt='%f')
np.savetxt("../SDECode/matlabPlots/ParameterEstimation/resultsTwo/ResultsFilterT_" + str(int(T)) + "Cova.txt", postCov, fmt='%f')

