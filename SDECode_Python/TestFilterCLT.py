from NISDE import *
from ParEst import *
from Hom import *
from multiprocessing import Pool
import functools
import matplotlib.pyplot as plt
from scipy.stats import norm 

sigma = 1.0
betaFilter = 1

# Set up the equation
eps = 0.05
sde = MSOrnUhl(eps)
sdeHomo = OrnUhl()
EM = EulerMaruyama(sde)
alpha = 1.0
IC = 0.0

# Set parameters
T = 50.
beta = 3.0
h = eps ** beta
N = int(round(T / h))
delta = 1.
print('T = ', T, 'h = ', h, 'N = ', round(T / h))

# Compute the homogenized coefficients
homParam = compute_homogeneous(sigma, alpha, 2.0 * np.pi, sde.p)
print(homParam)


def in_loop(seed):
    print seed
    # Generate solution
    y, dw = EM.solve(IC, N, T, alpha, sigma, seed=seed)
    filt, kt = filter_trajectory(y, delta, T, beta=betaFilter, type=1)
    par_est_filt = ParEst(filt, sdeHomo.grad_v_vect, h, is_vect=True)
    a_fi = par_est_filt.drift_mixed(data=y)
    return a_fi

nExp = 1000
p = Pool(min(nExp, 10))
aFilt = np.asarray(p.map(functools.partial(in_loop), range(nExp)))
p.close()

TLong = 5*T
NLong = int(round(TLong / h))
y, dw = EM.solve(IC, NLong, TLong, alpha, sigma, seed=nExp+1)
filt, kt = filter_trajectory(y, delta, TLong, beta=betaFilter, type=1)
parEstFilt = ParEst(filt, sdeHomo.grad_v_vect, h, is_vect=True)
aFi, num, den = parEstFilt.drift_mixed(data=y, return_exp=True)
den = den ** 2
varTh = 2 * sigma * num / den

plt.hist(np.sqrt(T) * (aFilt - homParam[0]), normed=True, bins=50)
empStdDev = np.sqrt(T) * aFilt.std()
x = np.linspace(norm.ppf(0.005, scale=empStdDev), norm.ppf(0.995, scale=empStdDev), 100)
plt.plot(x, norm.pdf(x, scale=empStdDev))
plt.show()

np.savetxt("../SDECode/matlabPlots/ParameterEstimation/resultsCLT/OUHom" + str(int(T)) + ".txt", [homParam[0]], fmt='%f')
np.savetxt("../SDECode/matlabPlots/ParameterEstimation/resultsCLT/OU" + str(int(T)) + ".txt", aFilt, fmt='%f')
np.savetxt("../SDECode/matlabPlots/ParameterEstimation/resultsCLT/OUVar" + str(int(T)) + ".txt", [varTh], fmt='%f')
