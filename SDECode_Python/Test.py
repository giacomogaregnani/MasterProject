from NISDE import *
from ParEst import *
from Hom import *
import matplotlib.pyplot as plt
import sys


# Set up the equation
eps = 0.05
sde = MSSixth(eps)
sdeHomo = Sixth()
EM = EulerMaruyama(sde)
alpha = 1.0
sigma = 1.0
param = np.array([alpha, sigma])
IC = 0

# Compute the homogenized coefficients
hom_param = compute_homogeneous(param, 2.0 * np.pi, sde.p)
print(hom_param)

# Set parameters
gamma = 1.3
beta = 2.0
zeta = beta - 0.5  # For zeta in (0, 1) in PaS07, fix here zeta in (beta-1, beta)
zetaFilter = 0.5
T = round(eps ** (-1.0 * gamma))
h = eps ** beta
N = int(round(T/h))
print('T = ', T, 'h = ', h)
deltaSub = int(eps ** (-1.0 * zeta))
deltaFilter = eps ** zetaFilter

tVec = np.arange(N+1) * h
nSub = N / deltaSub
tVecSub = np.arange(nSub+1) * h * deltaSub

# Initialize results
nExp = 100
A_filt = np.zeros(nExp)
A_sub = np.zeros(nExp)
A_interp = np.zeros(nExp)
A_bayes = np.zeros(nExp)
A_full = np.zeros(nExp)
A_mov = np.zeros(nExp)
var_bayes = np.zeros(nExp)

S_filt = np.zeros(nExp)
S_sub = np.zeros(nExp)
S_interp = np.zeros(nExp)
S_bayes = np.zeros(nExp)
S_full = np.zeros(nExp)
S_mov = np.zeros(nExp)

for i in range(0, nExp):
    print(i)

    # Generate solution
    print 'generating data...',
    Y, dw = EM.solve(IC, N, T, param)
    print(' data generated')

    # With subsampling
    print('sub'),
    sub = Y[::deltaSub]
    parEstSub = ParEst(sub, sdeHomo.grad_v_vect, deltaSub * h, is_vect=True)
    # A_sub[i] = parEstSub.drift()
    A_sub[i], S_sub[i] = parEstSub.drift_alternative(sdeHomo.lapl_v_vect, hom_param[1], estimate=True)
    print(A_sub[i], S_sub[i])

    # With subsampling and interpolation
    print('int'),
    tVecInterp = tVec[0:(np.size(sub)-1)*deltaSub]
    subPiecewise = np.interp(tVecInterp, tVecSub, sub)
    parEstSubPW = ParEst(subPiecewise, sdeHomo.grad_v_vect, h, is_vect=True)
    # A_interp[i] = parEstSubPW.drift(strat=True, sigma=param[1], lapl=sdeHomo.lapl_v_vect, Sigma=hom_param[1])
    A_interp[i], S_interp[i] = parEstSubPW.drift_alternative(sdeHomo.lapl_v_vect, hom_param[1], estimate=True, delta=1.0/(h*deltaSub))
    print(A_interp[i], S_interp[i])

    # With filtering
    filter, kt = filter_trajectory(Y, deltaFilter, T, beta=1)
    print('fil'),
    parEstFilt = ParEst(filter, sdeHomo.grad_v_vect, h, is_vect=True)
    # A_filt[i] = parEstFilt.drift(strat=True, sigma=param[1], lapl=sdeHomo.lapl_v_vect, Sigma=hom_param[1]) / (kt)
    # S_filt[i] = parEstFilt.diffusion() / (kt * h)
    A_filt[i], S_filt[i] = parEstFilt.drift_alternative(sdeHomo.lapl_v_vect, hom_param[1], estimate=True, delta=kt*h/2)
    print(A_filt[i], S_filt[i])

    # With a moving average
    movAv = moving_average(Y, int(round(deltaFilter / h)))
    print('mov'),
    parEstMov = ParEst(movAv, sdeHomo.grad_v_vect, h, is_vect=True)
    # A_mov[i] = parEstMov.drift(strat=True, sigma=param[1], lapl=sdeHomo.lapl_v_vect, Sigma=hom_param[1])
    A_mov[i], S_mov[i] = parEstMov.drift_alternative(sdeHomo.lapl_v_vect, hom_param[1], estimate=True, delta=1.0/deltaSub)
    print(A_mov[i], S_mov[i])

    # With nothing
    print('ful'),
    parEstFull = ParEst(Y, sdeHomo.grad_v_vect, h, is_vect=True)
    # A_full[i] = parEstFull.drift()
    A_full[i], S_full[i] = parEstFull.drift_alternative(sdeHomo.lapl_v_vect, hom_param[1], estimate=True)
    print(A_full[i], S_full[i])

    # Posterior
    A_bayes[i], var_bayes[i] = parEstFilt.driftBayesian(0.0, 1.0)

    if i == 0:
        plt.plot(tVec, Y)
        plt.plot(tVec, filter)
        plt.plot(tVec, movAv)
        plt.plot(tVecInterp, subPiecewise)
        plt.plot(tVecSub, Y[::deltaSub], '.')
        plt.legend(['obs', 'filter', 'movav', 'interp', 'subsampling'])
        plt.show()

print(['hom: ', hom_param[0], hom_param[1]])
print(['fil: ', 'A ', A_filt.mean(), A_filt.std(), 'S ', S_filt.mean(), S_filt.std()])
print(['mov: ', 'A ', A_mov.mean(), A_mov.std(), 'S ', S_mov.mean(), S_mov.std()])
print(['int: ', 'A ', A_interp.mean(), A_interp.std(), 'S ', S_interp.mean(), S_interp.std()])
print(['sub: ', 'A ', A_sub.mean(), A_sub.std(), 'S ', S_sub.mean(), S_sub.std()])
print(['ful: ', 'A ', A_full.mean(), A_full.std(), 'S ', S_full.mean(), S_full.std()])

if nExp == 1:
    sys.exit()

A_filt_sorted, A_filt_density = get_density(A_filt)
A_interp_sorted, A_interp_density = get_density(A_interp)
A_sub_sorted, A_sub_density = get_density(A_sub)
A_full_sorted, A_full_density = get_density(A_full)
A_mov_sorted, A_mov_density = get_density(A_mov)

plt.hist(A_filt, bins=30, normed=True, alpha=0.5)
plt.hist(A_sub, bins=30, normed=True, alpha=0.5)
plt.hist(A_full, bins=30, normed=True, alpha=0.5)
#plt.hist(A_interp, bins=20, normed=True, alpha=0.5)
#plt.hist(A_mov, bins=20, normed=True, alpha=0.5)
plt.plot(A_filt_sorted, A_filt_density)
plt.plot(A_sub_sorted, A_sub_density)
plt.plot(A_full_sorted, A_full_density)
#plt.plot(A_interp_sorted, A_interp_density)
#plt.plot(A_mov_sorted, A_mov_density)
plt.axvline(x=hom_param[0])
plt.axvline(x=param[0])
plt.legend(['filt', 'subsam', 'full']) #, 'interp', 'movav'])
plt.show()

S_filt_sorted, S_filt_density = get_density(S_filt)
S_interp_sorted, S_interp_density = get_density(S_interp)
S_sub_sorted, S_sub_density = get_density(S_sub)
S_full_sorted, S_full_density = get_density(S_full)
S_mov_sorted, S_mov_density = get_density(S_mov)

plt.hist(S_filt, bins=30, normed=True, alpha=0.5)
plt.hist(S_sub, bins=30, normed=True, alpha=0.5)
plt.hist(S_full, bins=30, normed=True, alpha=0.5)
#plt.hist(S_interp, bins=20, normed=True, alpha=0.5)
#plt.hist(S_mov, bins=20, normed=True, alpha=0.5)
plt.plot(S_filt_sorted, S_filt_density)
plt.plot(S_sub_sorted, S_sub_density)
plt.plot(S_full_sorted, S_full_density)
#plt.plot(S_interp_sorted, S_interp_density)
#plt.plot(S_mov_sorted, S_mov_density)
plt.axvline(x=hom_param[1])
plt.axvline(x=param[1])
plt.legend(['filt', 'subsam', 'full']) #, 'interp', 'movav'])
plt.show()

# Open a file for plot
np.savetxt("../SDECode/matlabPlots/ParameterEstimation/ResultsFilterHom.txt", hom_param, fmt='%f')
np.savetxt("../SDECode/matlabPlots/ParameterEstimation/ResultsFilter.txt", A_filt, fmt='%f')
np.savetxt("../SDECode/matlabPlots/ParameterEstimation/ResultsFilterSubs.txt", A_sub, fmt='%f')
np.savetxt("../SDECode/matlabPlots/ParameterEstimation/ResultsFilterFull.txt", A_full, fmt='%f')

np.savetxt("../SDECode/matlabPlots/ParameterEstimation/ResultsFilterDiff.txt", S_filt, fmt='%f')
np.savetxt("../SDECode/matlabPlots/ParameterEstimation/ResultsFilterSubsDiff.txt", S_sub, fmt='%f')
np.savetxt("../SDECode/matlabPlots/ParameterEstimation/ResultsFilterFullDiff.txt", S_full, fmt='%f')