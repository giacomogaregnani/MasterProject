from NISDE import *
from multiprocessing import Pool
from ParEst import *
import functools

sigma = 0.5
alpha = 3.
param = np.array([alpha, sigma])
IC = 0.0

# Set up the equation
eps = 0.05
sde = Sixth()
EM = EulerMaruyama(sde)

# Set parameters
betaFilter = 1
T = 50.
h = 0.05  # eps ** 2.5
N = int(round(T / h))
delta = 1


def generate_sol(seed):
    print seed
    y, dw = EM.solve(IC, N, T, param, seed=seed)
    filt, kt = filter_trajectory(y, delta, T, beta=betaFilter, type=1)
    return y[-1], filt[-1]

nExp = 11000
p = Pool(min(nExp, 11))
results = p.map(functools.partial(generate_sol), range(nExp))
yFin, zFin = np.asarray(zip(*results))
p.close()

# Write data to file
np.savetxt("../SDECode/matlabPlots/ParameterEstimation/GaussianTest_SI_HOM.txt", np.append(yFin, zFin), fmt='%f')
