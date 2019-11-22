from NISDE import *
from ParEst import *
from Hom import *
import matplotlib.pyplot as plt
import sys


# Set up the equation
eps = 0.1
sde = MSQuartic(eps)
sdeHomo = Quartic()
EM = EulerMaruyama(sde)
alpha = 1.0
sigma = .5
param = np.array([alpha, sigma])
IC = 0

# Compute the homogenized coefficients
hom_param = compute_homogeneous(param, 2.0 * np.pi, sde.p)
print(hom_param)

# Set parameters
betaFilter = 4.0
gamma = 1.5
beta = 2.0
zeta = beta - 0.33  # For zeta in (0, 1) in PaS07, fix here zeta in (beta-1, beta)
zetaFilter = 1.0
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
nExp = 200
A_filt = np.zeros(nExp)
A_sub = np.zeros(nExp)
A_full = np.zeros(nExp)

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
    A_sub[i] = parEstSub.drift()
    print(A_sub[i])

    # With filtering
    filter, kt = filter_trajectory(Y, deltaFilter, T, beta=betaFilter)
    print('fil'),
    parEstFilt = ParEst(filter, sdeHomo.grad_v_vect, h, is_vect=True)
    S_filt = parEstFilt.diffusion() * (2.0**(1.0/betaFilter))/(h*kt)
    A_filt[i] = parEstFilt.drift(strat=True, lapl=sdeHomo.lapl_v_vect, Sigma=S_filt)
    print(A_filt[i], S_filt)

    # With nothing
    print('ful'),
    parEstFull = ParEst(Y, sdeHomo.grad_v_vect, h, is_vect=True)
    A_full[i] = parEstFull.drift()
    print(A_full[i])

    # Posterior
    if i == 0:
        plt.plot(tVec, Y)
        plt.plot(tVec, filter)
        plt.plot(tVecSub, Y[::deltaSub], '.')
        plt.legend(['obs', 'filter', 'subsampling'])
        plt.show()

print(['hom: ', hom_param[0], hom_param[1]])
print(['fil: ', 'A ', A_filt.mean(), A_filt.std(), 'S '])
print(['sub: ', 'A ', A_sub.mean(), A_sub.std(), 'S '])
print(['ful: ', 'A ', A_full.mean(), A_full.std(), 'S '])

if nExp == 1:
    sys.exit()

A_filt_sorted, A_filt_density = get_density(A_filt)
A_sub_sorted, A_sub_density = get_density(A_sub)
A_full_sorted, A_full_density = get_density(A_full)

plt.hist(A_filt, bins=30, normed=True, alpha=0.5)
plt.hist(A_sub, bins=30, normed=True, alpha=0.5)
plt.hist(A_full, bins=30, normed=True, alpha=0.5)
plt.plot(A_filt_sorted, A_filt_density)
plt.plot(A_sub_sorted, A_sub_density)
plt.plot(A_full_sorted, A_full_density)
plt.axvline(x=hom_param[0])
plt.axvline(x=param[0])
plt.legend(['filt', 'subsam', 'full'])
plt.show()

# Open a file for plot
np.savetxt("../SDECode/matlabPlots/ParameterEstimation/ResultsFilterHomStrato.txt", hom_param, fmt='%f')
np.savetxt("../SDECode/matlabPlots/ParameterEstimation/ResultsFilterStrato.txt", A_filt, fmt='%f')
np.savetxt("../SDECode/matlabPlots/ParameterEstimation/ResultsFilterSubsStrato.txt", A_sub, fmt='%f')
np.savetxt("../SDECode/matlabPlots/ParameterEstimation/ResultsFilterFullStrato.txt", A_full, fmt='%f')