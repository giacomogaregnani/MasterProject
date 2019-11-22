from NISDE import *
from ParEst import *
from Hom import *
import matplotlib.pyplot as plt
from multiprocessing import Pool

# Set up the equation
eps = 0.1
sde = MSQuartic(eps)
sdeHomo = Quartic()
EM = EulerMaruyama(sde)
alpha = 1.0
sigma = .6
param = np.array([alpha, sigma])
IC = 0.0

# Compute the homogenized coefficients
hom_param = compute_homogeneous(param, 2.0 * np.pi, sde.p)
print(hom_param)

# Set parameters
betaFilt = 1.0
gamma = 2.0
beta = 3.0
zetaFilt = 1.0
T = round(eps ** (-1.0 * gamma))
h = eps ** beta
N = int(round(T/h))
print('T = ', T, 'h = ', h)
deltaFilt = eps ** zetaFilt

nExp = 500

def inLoop(seed):
    print seed
    Y, dW = EM.solve(IC, N, T, param, seed=seed)
    filter, kt = filter_trajectory(Y, deltaFilt, T, beta=betaFilt)
    ftg = 1.0 / (N + 1) * np.sum(sdeHomo.grad_v_vect(Y) ** 2)
    ffg = 1.0 / (N + 1) * np.sum(sdeHomo.grad_v_vect(filter) ** 2)
    ftl = 1.0 / (N + 1) * np.sum(sdeHomo.lapl_v_vect(Y))
    ffl = 1.0 / (N + 1) * np.sum(sdeHomo.lapl_v_vect(filter))
    return ftg, ffg, ftl, ffl

p = Pool()
finalTrueGrad, finalFiltGrad, finalTrueLapl, finalFiltLapl = np.asarray(zip(*p.map(inLoop, range(nExp))))

print(finalTrueGrad.mean(), sigma/alpha * finalTrueLapl.mean())

finalTrueGradSorted, finalTrueGradDensity = get_density(finalTrueGrad)
finalFiltGradSorted, finalFiltGradDensity = get_density(finalFiltGrad)

plt.hist(finalTrueGrad, bins=30, density=True, alpha=0.5)
plt.hist(finalFiltGrad, bins=30, density=True, alpha=0.5)
plt.plot(finalTrueGradSorted, finalTrueGradDensity)
plt.plot(finalFiltGradSorted, finalFiltGradDensity)
plt.legend(['true', 'filter'])
plt.show()

plt.plot(finalTrueGradSorted)
plt.plot(finalFiltGradSorted)
plt.show()

finalTrueLaplSorted, finalTrueLaplDensity = get_density(finalTrueLapl)
finalFiltLaplSorted, finalFiltLaplDensity = get_density(finalFiltLapl)

plt.hist(finalTrueLapl, bins=30, density=True, alpha=0.5)
plt.hist(finalFiltLapl, bins=30, density=True, alpha=0.5)
plt.plot(finalTrueLaplSorted, finalTrueLaplDensity)
plt.plot(finalFiltLaplSorted, finalFiltLaplDensity)
plt.legend(['true', 'filter'])
plt.show()

plt.plot(finalTrueLaplSorted)
plt.plot(finalFiltLaplSorted)
plt.show()

