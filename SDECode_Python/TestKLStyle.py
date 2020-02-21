from NISDE import *
from ParEst import *
from Hom import *
import matplotlib.pyplot as plt
from multiprocessing import Pool
import functools

sigmaVec = [0.7]
betaFilterVec = [1]

nParam = 6
alpha = np.zeros(nParam)
alpha[0] = -0.5
alpha[1] = -0.25
alpha[2] = -0.125
alpha[3] = 0.125
alpha[4] = 0.25
alpha[5] = 0.5

# Set up the equation
eps = 0.1
sde = MSCheb(eps, nParam)
sdeHomo = Cheb(nParam)
EM = EulerMaruyama(sde)

np.savetxt("../SDECode/matlabPlots/ParameterEstimation/resultsKL/ResultsAlpha.txt", alpha, fmt='%f')

for sigma in sigmaVec:

    IC = 0.0

    # Compute the homogenized coefficients
    hom_param = compute_homogeneous(sigma, alpha, 2.0 * np.pi, sde.p)
    print(hom_param)

    np.savetxt("../SDECode/matlabPlots/ParameterEstimation/resultsKL/ResultsHom" + "_s" + str(sigma) + ".txt",
               hom_param, fmt='%f')

    for betaFilter in betaFilterVec:

        # Set parameters
        gammaTime = -np.log(100000) / np.log(eps)
        beta = 3.0
        T = round(eps ** (-1.0 * gammaTime))
        h = eps ** beta
        print('T = ', T, 'h = ', h, 'N = ', round(T/h))

        nZeta = 1
        zetaSubsVec = [0.66]
        zetaFilterVec = [0.0]

        A_filt_means = np.zeros((nZeta, nParam+1))
        A_sub_means = np.zeros((nZeta, nParam+1))
        A_noth_means = np.zeros((nZeta, nParam))

        for i in range(0, nZeta):

            N = int(round(T/h))

            delta = eps ** zetaSubsVec[i]
            nSub = int(delta / h)
            tVec = np.arange(N+1) * h
            deltaFilter = eps ** zetaFilterVec[i]

            # Initialize results
            nExp = 1

            # plot = False
            # for seed in range(0, nExp):
            def in_loop(seed, plot=False):
                print seed
                y, dw = EM.solve(IC, N, T, alpha, sigma, seed=seed)
                par_est_nothing = ParEst(y, sdeHomo.grad_v, h)
                sub = y[::nSub]
                par_est_sub = ParEst(sub, sdeHomo.grad_v, delta, is_vect=True)
                filt, kt = filter_trajectory(y, deltaFilter, T, beta=betaFilter, type=1)
                par_est_filt = ParEst(filt, sdeHomo.grad_v, h, is_vect=True)

                a_nothing = par_est_nothing.drift_multi(nParam)
                a_fi = par_est_filt.drift_multi_mixed(nParam, data=y, tilde=False)
                a_sub = par_est_sub.drift_multi(nParam)

                if plot and seed == 0:
                    plt.plot(y)
                    plt.plot(filt)
                    plt.show()

                return a_sub, a_fi, a_nothing

            p = Pool(min(nExp, 1))
            results = p.map(functools.partial(in_loop, plot=False), range(nExp))
            aSub, aFilt, aNoth = np.asarray(zip(*results))
            p.close()

            # Compute means
            aFiltMean = np.zeros(nParam)
            aSubMean = np.zeros(nParam)
            aNothMean = np.zeros(nParam)
            for j in range(0, nExp):
                aFiltMean += aFilt[j]
                aSubMean += aSub[j]
                aNothMean += aNoth[j]
            aFiltMean /= nExp
            aSubMean /= nExp
            aNothMean /= nExp

            print(['hom: ', 'A ', hom_param[0:-1]])
            print(['fil: ', 'A ', aFiltMean])
            print(['not: ', 'A ', aNothMean])
            print(['sub: ', 'A ', aSubMean])

            A_filt_means[i] = np.append(aFiltMean, [deltaFilter])
            A_sub_means[i] = np.append(aSubMean, [delta])
            A_noth_means[i] = aNothMean

        # Write data to file
        np.savetxt("../SDECode/matlabPlots/ParameterEstimation/resultsKL/ResultsFilter" + "_s" + str(sigma) + "_b" + str(betaFilter) + ".txt", A_filt_means, fmt='%f')
    np.savetxt("../SDECode/matlabPlots/ParameterEstimation/resultsKL/ResultsSub" + "_s" + str(sigma) + ".txt" , A_sub_means, fmt='%f')
    np.savetxt("../SDECode/matlabPlots/ParameterEstimation/resultsKL/ResultsNot" + "_s" + str(sigma) + ".txt" , A_noth_means, fmt='%f')
