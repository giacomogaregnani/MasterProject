from NISDE import *
from ParEst import *
from Hom import *
from multiprocessing import Pool
import functools

sigmaVec = [.5, .7, 1.]
betaFilterVec = np.arange(10) + 1

# Set up the equation
eps = 0.1
sde = MSQuartic(eps)
sdeHomo = Quartic()
EM = EulerMaruyama(sde)
alpha = 1.0
IC = 0.0

# Set parameters
gammaTime = -np.log(1000) / np.log(eps)
beta = 3.0
T = round(eps ** (-1.0 * gammaTime))
h = eps ** beta
N = int(round(T / h))
delta = eps ** 0.5

print('T = ', T, 'h = ', h, 'N = ', round(T / h))

for sigma in sigmaVec:

    # Compute the homogenized coefficients
    hom_param = compute_homogeneous(sigma, alpha, 2.0 * np.pi, sde.p)
    print(hom_param)

    np.savetxt("../SDECode/matlabPlots/ParameterEstimation/resultsBeta/ResultsHom" + "_s" + str(sigma) + ".txt",
               [hom_param[0]], fmt='%f')

    aMeans = np.zeros(10)

    for i in range(0, 10):

        betaFilter = betaFilterVec[i]
        print('beta = ', betaFilter)

        # Initialize results
        nExp = 10

        def in_loop(seed):
            print seed
            # Generate solution
            y, dw = EM.solve(IC, N, T, alpha, sigma, seed=seed)
            filt, kt = filter_trajectory(y, delta, T, beta=betaFilter, type=1)
            par_est_filt = ParEst(filt, sdeHomo.grad_v_vect, h, is_vect=True)

            a_fi = par_est_filt.drift_mixed(data=y)

            return a_fi

        p = Pool(min(nExp, 10))
        aFilt = np.asarray(p.map(functools.partial(in_loop), range(nExp)))
        p.close()

        aMeans[i] = aFilt.mean()

        print(['hom: ', 'A ', hom_param[0]])
        print(['fil: ', 'A ', aFilt.mean()])

    # Write data to file
    np.savetxt("../SDECode/matlabPlots/ParameterEstimation/resultsBeta/ResultsFilter" + "_s" + str(sigma)
               + ".txt", aMeans, fmt='%f')
