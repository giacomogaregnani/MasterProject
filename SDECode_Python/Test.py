from NISDE import *
from ParEst import *
from Hom import *
import matplotlib.pyplot as plt
import sys
from multiprocessing import Pool
import functools


# Set up the equation
eps = 0.1
sde = MSOrnUhl(eps)
sdeHomo = OrnUhl()
EM = EulerMaruyama(sde)
alpha = 1.0
sigma = 1.0
param = np.array([alpha, sigma])
IC = 0.0

# Compute the homogenized coefficients
hom_param = compute_homogeneous(param, 2.0 * np.pi, sde.p)
print(hom_param)

# Set parameters
gamma = 2.0
betaFilter = 1.0
beta = 3.0

for i in range(0, 11):

    print i

    zetaSub = i*0.1
    zetaFilter = i*0.1
    T = round(eps ** (-1.0 * gamma))
    h = eps ** beta
    N = int(round(T/h))
    print('T = ', T, 'h = ', h)
    deltaSub = eps ** zetaSub
    nSub = int(round(deltaSub / h))
    deltaFilter = eps ** zetaFilter

    tVec = np.arange(N+1) * h

    # Initialize results
    nExp = 400

    def inLoop(seed, strato=False):
        # print(seed)

        # Generate solution
        Y, dw = EM.solve(IC, N, T, param, seed=seed)
        # With subsampling
        sub = Y[::nSub]
        parEstSub = ParEst(sub, sdeHomo.grad_v_vect, deltaSub, is_vect=True)
        # With filtering
        filter, kt = filter_trajectory(Y, deltaFilter, T, beta=betaFilter)
        parEstFilt = ParEst(filter, sdeHomo.grad_v_vect, h, is_vect=True)
        # With nothing
        parEstFull = ParEst(Y, sdeHomo.grad_v_vect, h, is_vect=True)
        if strato:
            Sfu = parEstFull.diffusion()
            Afu = parEstFull.drift()
            Sfi = parEstFilt.diffusion() / (kt*h/(2.0**(1.0/betaFilter)))
            # Afi = parEstFilt.drift() / (kt*h/(2.0**(1.0/betaFilter)))
            Afi = parEstFilt.drift(strat=True, Sigma=Sfi, lapl=sdeHomo.lapl_v_vect)
            Ss = parEstSub.diffusion()
            As = parEstSub.drift()
        else:
            Afi, Sfi = parEstFilt.drift_alternative(sdeHomo.lapl_v_vect, delta=kt*h/(2.0**(1.0/betaFilter)))
            Afu, Sfu = parEstFull.drift_alternative(sdeHomo.lapl_v_vect)
            As, Ss = parEstSub.drift_alternative(sdeHomo.lapl_v_vect)
        return As, Ss, Afi, Sfi, Afu, Sfu


    p = Pool(11)
    results = p.map(functools.partial(inLoop, strato=False), range(nExp))
    A_sub, S_sub, A_filt, S_filt, A_full, S_full = np.asarray(zip(*results))

    print(['hom: ', 'A ', hom_param[0], 'S ', hom_param[1]])
    print(['fil: ', 'A ', A_filt.mean(), A_filt.std(), 'S ', S_filt.mean(), S_filt.std()])
    print(['sub: ', 'A ', A_sub.mean(), A_sub.std(), 'S ', S_sub.mean(), S_sub.std()])
    print(['ful: ', 'A ', A_full.mean(), A_full.std(), 'S ', S_full.mean(), S_full.std()])

if nExp == 1:
    sys.exit()

A_filt_sorted, A_filt_density = get_density(A_filt)
A_full_sorted, A_full_density = get_density(A_full)
A_sub_sorted, A_sub_density = get_density(A_sub)


plt.hist(A_filt, bins=30, density=True, alpha=0.5)
plt.hist(A_sub, bins=30, density=True, alpha=0.5)
# plt.hist(A_full, bins=30, density=True, alpha=0.5)
plt.plot(A_filt_sorted, A_filt_density)
plt.plot(A_sub_sorted, A_sub_density)
# plt.plot(A_full_sorted, A_full_density)
plt.axvline(x=hom_param[0])
plt.axvline(x=param[0])
plt.legend(['filt', 'subsam'])
plt.show()

S_filt_sorted, S_filt_density = get_density(S_filt)
S_sub_sorted, S_sub_density = get_density(S_sub)
S_full_sorted, S_full_density = get_density(S_full)

plt.hist(S_filt, bins=30, density=True, alpha=0.5)
plt.hist(S_sub, bins=30, density=True, alpha=0.5)
# plt.hist(S_full, bins=30, density=True, alpha=0.5)
plt.plot(S_filt_sorted, S_filt_density)
plt.plot(S_sub_sorted, S_sub_density)
# plt.plot(S_full_sorted, S_full_density)
plt.axvline(x=hom_param[1])
plt.axvline(x=param[1])
plt.legend(['filt', 'subsam']) #, 'full'])
plt.show()

# Write data to file
np.savetxt("../SDECode/matlabPlots/ParameterEstimation/ResultsFilterHom.txt", hom_param, fmt='%f')
np.savetxt("../SDECode/matlabPlots/ParameterEstimation/ResultsFilter.txt", A_filt, fmt='%f')
np.savetxt("../SDECode/matlabPlots/ParameterEstimation/ResultsFilterSubs.txt", A_sub, fmt='%f')
np.savetxt("../SDECode/matlabPlots/ParameterEstimation/ResultsFilterFull.txt", A_full, fmt='%f')

np.savetxt("../SDECode/matlabPlots/ParameterEstimation/ResultsFilterDiff.txt", S_filt, fmt='%f')
np.savetxt("../SDECode/matlabPlots/ParameterEstimation/ResultsFilterSubsDiff.txt", S_sub, fmt='%f')
np.savetxt("../SDECode/matlabPlots/ParameterEstimation/ResultsFilterFullDiff.txt", S_full, fmt='%f')