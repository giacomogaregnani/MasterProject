from NISDE import *
from ParEst import *
from Hom import *
import matplotlib.pyplot as plt

epsVec = np.array([0.1, 0.05, 0.025, 0.0125, 0.00625, 0.003125])
error = np.zeros(np.size(epsVec))

i = 0

for eps in epsVec:

    sde = MSOrnUhl(eps)
    EM = EulerMaruyama(sde)
    alpha = 1.0
    sigma = 1.0
    param = np.array([alpha, sigma])
    IC = 0.0

    # Set parameters
    betaFilt = 5.0
    gamma = 3.0
    beta = 2.0
    zetaFilt = 1.0
    T = 20.0
    h = np.min(epsVec) ** beta
    N = int(round(T/h))
    print('T = ', T, 'h = ', h)
    deltaFilt = eps ** zetaFilt

    Y, dW = EM.solve(IC, N, T, param, seed=1)
    filter, kt = filter_trajectory(Y, deltaFilt, T, beta=betaFilt)

    tVec = np.arange(N+1) * h
    # plt.plot(tVec, Y)
    # plt.plot(tVec, filter)
    # plt.show()

    error[i] = np.max(np.abs((filter - Y))) / T
    i += 1

plt.plot(epsVec, error, '-o')
plt.show()