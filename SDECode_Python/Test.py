from NISDE import *
from ParEst import *
from Hom import *
import matplotlib.pyplot as plt

# Set up the equation
eps = 0.04
eqn = MSOrnUhl(eps)
EM = EulerMaruyama(eqn)
param = np.array([1.0, 0.5])
IC = 0.0

# Compute the homogenized coefficients
hom_param = compute_homogeneous(param, 2*np.pi, eqn.grad_p)
print(hom_param)

# Set final time and moving average window
gamma = 3
beta = 2
zeta = 1.5
N = int(pow(eps, -1.0 * gamma))
h = pow(eps, beta)
T = h * N
delta = int(pow(eps, -zeta))

# Generate solution and compute its moving average and plot
Y = EM.solve(IC, N, T, param)
avg = moving_average(Y, delta)
plt.plot(Y)
plt.plot(avg)
plt.show()

# Estimate parameters
parEst = ParEst(Y, eqn.grad_v, h)
A = parEst.drift()
S = parEst.diffusion()
print([A, S])

# With subsampling
sub = Y[::delta]
parEstSub = ParEst(sub, eqn.grad_v, delta*h)
A_sub = parEstSub.drift()
S_sub = parEstSub.diffusion()
print([A_sub, S_sub])

# With averaging
parEstAvg = ParEst(avg, eqn.grad_v, h)
A_avg = parEstAvg.drift()
S_avg = delta * parEstAvg.diffusion()
print([A_avg, S_avg])
