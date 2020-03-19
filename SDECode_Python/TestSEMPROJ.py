from NISDE import *
from ParEst import *
from Hom import *
from multiprocessing import Pool
import functools
import matplotlib.pyplot as plt

sigma = 0.5

# Set up the equation
eps = 0.1
sde = MSQuartic(eps)
sdeHomo = Quartic()
EM = EulerMaruyama(sde)
alpha = 1.0
IC = 0.0

# Set parameters
T = 2000.
beta = 3.0
h = eps ** beta
N = int(round(T / h))
delta = eps ** 0.5
nRatio = int(delta / h)
print('T = ', T, 'h = ', h, 'N = ', round(T / h))

# Compute the homogenized coefficients
homParam = compute_homogeneous(sigma, alpha, 2.0 * np.pi, sde.p)
print(homParam)


def in_loop(seed):
    print seed
    # Generate solution
    y, dw = EM.solve(IC, N, T, alpha, sigma, seed=seed)
    y_sub = y[::nRatio]
    par_est_subs = ParEst(y_sub, sdeHomo.grad_v_vect, delta, is_vect=True)
    a_subs = par_est_subs.drift()
    plt.plot(y)
    plt.show()
    return a_subs

nExp = 1
p = Pool(min(nExp, 10))
aSubs = np.asarray(p.map(functools.partial(in_loop), range(nExp)))
p.close()

print aSubs.mean()
