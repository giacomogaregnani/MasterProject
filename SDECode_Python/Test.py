from NISDE import *
from ParEst import *
from Hom import *
import matplotlib.pyplot as plt
from multiprocessing import Pool
import functools
from scipy.stats import norm

sigmaVec = [.4]
betaVec = [1, 5, 10]

# Set up the equation
eps = 0.1
sde = MSQuartic(eps)
sdeHomo = Quartic()
EM = EulerMaruyama(sde)

for sigma in sigmaVec:

    alpha = 1.0
    param = np.array([alpha, sigma])
    IC = 0.0

    # Compute the homogenized coefficients
    hom_param = compute_homogeneous(param, 2.0 * np.pi, sde.p)
    print(hom_param)

    np.savetxt("../SDECode/matlabPlots/ParameterEstimation/resultsQU_epsBig/ResultsHom" + "_s" + str(sigma) + ".txt", hom_param, fmt='%f')

    for betaFilter in betaVec:

        # Set parameters
        gammaTime = -np.log(100) / np.log(eps)
        beta = 2.5
        T = round(eps ** (-1.0 * gammaTime))
        h = eps ** beta
        print('T = ', T, 'h = ', h, 'N = ', round(T/h))

        nZeta = 11
        zMax = 1.0
        zIncrement = zMax / (nZeta - 1)

        A_filt_means = np.zeros(nZeta)
        A_filt_stds = np.zeros(nZeta)
        A_sub_means = np.zeros(nZeta)
        A_sub_stds = np.zeros(nZeta)
        S_filt_means = np.zeros(nZeta)
        S_filt_stds = np.zeros(nZeta)
        S_sub_means = np.zeros(nZeta)
        S_sub_stds = np.zeros(nZeta)

        supnSub = int(eps**(1.0-beta))+1
        infnSub = int(eps**(-beta))
        nSubRange = np.logspace(-1.0*beta, 1.0-beta, num=nZeta, base=eps)

        for i in range(0, nZeta):

            N = int(round(T/h))

            zeta = beta + np.log(nSubRange[i]) / np.log(eps)

            print(betaFilter, zeta)

            delta = eps ** zeta
            nSub = int(nSubRange[i])
            tVec = np.arange(N+1) * h

            # Initialize results
            nExp = 55

            def in_loop(seed, tilde=False, plot=False):
                print seed
                # Generate solution
                y, dw = EM.solve(IC, N, T, param, seed=seed)
                plt.show()
                sub = y[::nSub]
                par_est_sub = ParEst(sub, sdeHomo.grad_v_vect, delta, is_vect=True)
                filt, kt = filter_trajectory(y, delta, T, beta=betaFilter, type=1)
                par_est_filt = ParEst(filt, sdeHomo.grad_v_vect, h, is_vect=True)
                if tilde:
                    a_fi, s_fi = par_est_filt.drift_tilde(sdeHomo.lapl_v_vect, delta=kt * h / (2.0 ** (1.0 / betaFilter)))
                    a_sub, s_sub = par_est_sub.drift_tilde(sdeHomo.lapl_v_vect)
                else:
                    s_fi = par_est_filt.diffusion() / (kt*h / (2.0**(1.0/betaFilter)))
                    a_fi = par_est_filt.drift_other_data(data=y)
                    # verif1, verif2 = par_est_filt.formula_verification(data=y, grad_p=sde.grad_p,
                    #                                                    eps=eps, dw=dw, sigma=sigma)
                    # print(verif1, verif2, hom_param[0]-alpha)
                    s_sub = par_est_sub.diffusion()
                    a_sub = par_est_sub.drift()

                if plot and seed == 0:
                    plt.plot(y)
                    plt.plot(filt)
                    plt.show()

                return a_sub, s_sub, a_fi, s_fi  # , verif1, verif2

            p = Pool(min(nExp, 11))
            results = p.map(functools.partial(in_loop, tilde=False, plot=False), range(nExp))
            # aSub, sSub, aFilt, sFilt, verOne, verTwo = np.asarray(zip(*results))
            aSub, sSub, aFilt, sFilt = np.asarray(zip(*results))

            p.close()

            # res = np.sort(verTwo)
            # plt.hist(np.sqrt(T) * res, density=True, alpha=0.5, bins=20)
            # plt.axvline(x=0)
            # plt.plot(res, norm.pdf(res, 0, res.std()))
            # ax = plt.gca()
            # ax.set(xlim=(-8, 8))
            # plt.show()

            # res = np.sort(np.sqrt(T) * (aFilt - hom_param[0]))
            # plt.hist(res, density=True, alpha=0.5, bins=20)
            # # plt.axvline(x=hom_param[0])
            # plt.plot(res, norm.pdf(res, 0, res.std()))
            # plt.plot(res, norm.pdf(res, 0, np.sqrt(2*hom_param[0])))
            # plt.show()

            # res = np.sort(np.sqrt(T) * (aSub - hom_param[0]))
            # plt.hist(res, density=True, alpha=0.5, bins=20)
            # # plt.axvline(x=hom_param[0])
            # plt.plot(res, norm.pdf(res, 0, res.std()))
            # plt.plot(res, norm.pdf(res, 0, np.sqrt(2*hom_param[0])))
            # plt.show()

            # res = np.sort(np.sqrt(T) * (sFilt - hom_param[1]))
            # plt.hist(res, density=True, alpha=0.5, bins=20)
            # # plt.axvline(x=hom_param[1])
            # plt.plot(res, norm.pdf(res, 0, res.std()))
            # plt.show()

            # res = np.sort(np.sqrt(T) * (sSub - hom_param[1]))
            # plt.hist(res, density=True, alpha=0.5, bins=20)
            # # plt.axvline(x=hom_param[1])
            # plt.plot(res, norm.pdf(res, 0, res.std()))
            # plt.show()

            A_filt_means[i] = aFilt.mean()
            A_filt_stds[i] = aFilt.std()
            A_sub_means[i] = aSub.mean()
            A_sub_stds[i] = aSub.std()
            S_filt_means[i] = sFilt.mean()
            S_filt_stds[i] = sFilt.std()
            S_sub_means[i] = sSub.mean()
            S_sub_stds[i] = sSub.std()

            print(['hom: ', 'A ', hom_param[0], 'S ', hom_param[1]])
            print(['fil: ', 'A ', aFilt.mean(), aFilt.std(), 'S ', sFilt.mean(), sFilt.std()])
            print(['sub: ', 'A ', aSub.mean(), aSub.std(), 'S ', sSub.mean(), sSub.std()])

        # Write data to file
        np.savetxt("../SDECode/matlabPlots/ParameterEstimation/resultsQU_epsBig/ResultsFilter" + "_s" + str(sigma) + "_b" + str(betaFilter) + ".txt", np.append(A_filt_means, A_filt_stds), fmt='%f')
    np.savetxt("../SDECode/matlabPlots/ParameterEstimation/resultsQU_epsBig/ResultsSub" + "_s" + str(sigma) + ".txt" , np.append(A_sub_means, A_sub_stds), fmt='%f')
