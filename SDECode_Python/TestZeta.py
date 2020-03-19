from NISDE import *
from ParEst import *
from Hom import *
import matplotlib.pyplot as plt
from multiprocessing import Pool
import functools


sigmaVec = [.5]  # , .7, 1.]
betaFilterVec = [1]

# Set up the equation
eps = 0.1
sde = MSOrnUhl(eps)
sdeHomo = OrnUhl()
EM = EulerMaruyama(sde)
alpha = 1.0
IC = 0.0

for sigma in sigmaVec:

    # Compute the homogenized coefficients
    hom_param = compute_homogeneous(sigma, alpha, 2.0 * np.pi, sde.p)
    print(hom_param)

    np.savetxt("../SDECode/matlabPlots/ParameterEstimation/resultsOU_OtherFilter/ResultsHom" + "_s" + str(sigma) + ".txt",
               [hom_param[0]], fmt='%f')

    for betaFilter in betaFilterVec:

        # Set parameters
        gammaTime = -np.log(10000) / np.log(eps)
        beta = 3.0
        T = round(eps ** (-1.0 * gammaTime))
        h = eps ** beta
        print('T = ', T, 'h = ', h, 'N = ', round(T/h))

        nZeta = 10
        zMax = 1.0
        zetaVec = np.zeros(nZeta)

        A_filt_means = np.zeros(nZeta)
        A_filt_stds = np.zeros(nZeta)
        A_sub_means = np.zeros(nZeta)
        A_sub_stds = np.zeros(nZeta)

        supnSub = int(eps**(1.0-beta))+1
        infnSub = int(eps**(-beta))
        nSubRange = np.logspace(-1.0*beta, 1.0-beta, num=nZeta, base=eps)

        for i in range(0, nZeta):

            N = int(round(T/h))

            zeta = beta + np.log(nSubRange[i]) / np.log(eps)
            delta = eps ** zeta
            nSub = int(nSubRange[i])
            tVec = np.arange(N+1) * h
            zetaVec[i] = zeta

            # Initialize results
            nExp = 10

            def in_loop(seed, plot=False):
                print seed
                # Generate solution
                y, dw = EM.solve(IC, N, T, alpha, sigma, seed=seed)
                sub = y[::nSub]
                par_est_sub = ParEst(sub, sdeHomo.grad_v_vect, delta, is_vect=True)
                filt, kt = filter_trajectory(y, delta, T, beta=betaFilter, type=1)
                par_est_filt = ParEst(filt, sdeHomo.grad_v_vect, h, is_vect=True)

                a_fi = par_est_filt.drift_mixed(data=y)
                a_sub = par_est_sub.drift()

                if plot and seed == 0:
                    plt.plot(y)
                    plt.plot(filt)
                    plt.show()

                return a_sub, a_fi

            p = Pool(min(nExp, 10))
            results = p.map(functools.partial(in_loop, plot=False), range(nExp))
            aSub, aFilt, = np.asarray(zip(*results))
            p.close()

            A_filt_means[i] = aFilt.mean()
            A_filt_stds[i] = aFilt.std()
            A_sub_means[i] = aSub.mean()
            A_sub_stds[i] = aSub.std()

            print(['hom: ', 'A ', hom_param[0]])
            print(['fil: ', 'A ', aFilt.mean(), aFilt.std()])
            print(['sub: ', 'A ', aSub.mean(), aSub.std()])

        # Write data to file
        np.savetxt("../SDECode/matlabPlots/ParameterEstimation/resultsOU_OtherFilter/ResultsFilter" + "_s" + str(sigma) + "_b" + str(betaFilter) + ".txt",
                   np.append(A_filt_means, A_filt_stds), fmt='%f')
    np.savetxt("../SDECode/matlabPlots/ParameterEstimation/resultsOU_OtherFilter/ResultsSub" + "_s" + str(sigma) + ".txt",
               np.append(A_sub_means, A_sub_stds), fmt='%f')
    np.savetxt("../SDECode/matlabPlots/ParameterEstimation/resultsOU_OtherFilter/zetas.txt", zetaVec, fmt='%f')
