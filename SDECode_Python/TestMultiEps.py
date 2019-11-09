from NISDE import *
from ParEst import *
from Hom import *
import matplotlib.pyplot as plt


epsVec = np.array([0.2, 0.18, 0.16, 0.14, 0.12, 0.1, 0.08, 0.06, 0.04])
A_filt_mean = np.zeros(np.size(epsVec))
A_sub_mean = np.zeros(np.size(epsVec))
A_full_mean = np.zeros(np.size(epsVec))

for j in range(0, np.size(epsVec)):
    # Set up the equation
    eps = epsVec[j]
    sde = MSQuartic(eps)
    sdeHomo = Quartic()
    EM = EulerMaruyama(sde)
    param = np.array([1.0, 0.5])
    IC = 0.0

    # Compute the homogenized coefficients
    hom_param = compute_homogeneous(param, 2 * np.pi, sde.p)
    print(hom_param)

    # Set final time and moving average window
    gamma = 3.5
    beta = 2.0
    zeta = beta # For zeta in (0, 1) in PaS07, fix here zeta in (beta-1, beta)
    zetaFilter = 0.7
    T = round(eps ** (-1.0 * gamma))
    h = eps ** beta
    N = int(round(T/h))
    print(T, h, N)
    deltaSub = int(eps ** (-1.0 * zeta))
    deltaFilter = eps ** zetaFilter

    nExp = 1
    A_filt = np.zeros(nExp)
    A_sub = np.zeros(nExp)
    A_bayes = np.zeros(nExp)
    A_full = np.zeros(nExp)
    var_bayes = np.zeros(nExp)

    tVec = np.arange(N + 1) * h
    nSub = N / deltaSub
    tVecSub = np.arange(nSub+1) * h * deltaSub

    for i in range(0, nExp):
        print(i)

        # Generate solution and compute its filtered version
        Y = EM.solve(IC, N, T, param)
        avg, kt = filter_trajectory(Y, deltaFilter, T)

        # With subsampling
        sub = Y[::deltaSub]
        parEstSub = ParEst(sub, sdeHomo.grad_v_vect, deltaSub * h, is_vect=True)
        A_sub[i] = parEstSub.drift()

        # With filtering
        parEstFilt = ParEst(avg, sdeHomo.grad_v_vect, h, is_vect=True)
        A_filt[i] = parEstFilt.drift()

        # With nothing
        parEstFull = ParEst(Y, sdeHomo.grad_v_vect, h, is_vect=True)
        A_full[i] = parEstFull.drift()

        # Posterior
        A_bayes[i], var_bayes[i] = parEstFilt.driftBayesian(0.0, 1.0)
        # A_bayes[i] *= kt

    A_filt_mean[j] = A_filt.mean()
    A_sub_mean[j] = A_sub.mean()
    A_full_mean[j] = A_full.mean()

plt.loglog(epsVec, A_filt_mean)
plt.loglog(epsVec, A_sub_mean)
plt.loglog(epsVec, A_full_mean)
plt.loglog(epsVec, epsVec)
plt.legend(['filt', 'sub', 'full'])
plt.show()
