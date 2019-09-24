from NISDE import *
from ParEst import *
from Hom import *
import matplotlib.pyplot as plt

# Set up the equation
eps = 0.05
eqn = MSOrnUhl(eps)
eqn_homo = OrnUhl()
EM = EulerMaruyama(eqn)
param = np.array([1.0, 0.5])
IC = 0.0

# Compute the homogenized coefficients
hom_param = compute_homogeneous(param, 2*np.pi, eqn.grad_p)

# Set final time and moving average window
gamma = 4.0
beta = 2.0
zeta = 1.1
N = int(pow(eps, -1.0 * gamma))
h = pow(eps, beta)
T = h * N
delta = int(pow(eps, -1.0 * zeta))


n_exp = 500
A_avg = np.zeros(n_exp)
S_avg = np.zeros(n_exp)
A_sub = np.zeros(n_exp)
S_sub = np.zeros(n_exp)

for i in range(0, n_exp):
    print(i)

    # Generate solution and compute its moving average
    Y = EM.solve(IC, N, T, param)
    avg = moving_average(Y, delta)

    if i == n_exp-1:
        plt.plot(Y)
        plt.plot(avg)
        plt.show()

    # With subsampling
    # print('subsampling')
    sub = Y[::delta]
    parEstSub = ParEst(sub, eqn_homo.grad_v_vect, delta*h, is_vect=True)
    A_sub[i] = parEstSub.drift()
    S_sub[i] = parEstSub.diffusion()

    # With averaging
    # print('averaging')
    parEstAvg = ParEst(avg, eqn_homo.grad_v_vect, h, is_vect=True)
    A_avg[i] = delta * parEstAvg.drift()
    S_avg[i] = delta * parEstAvg.diffusion()
    # print([A_avg, S_avg])

print('results:')
print(hom_param)

print([A_avg.mean(), S_avg.mean()])
print([A_sub.mean(), S_sub.mean()])

A_avg_sorted, A_avg_density = get_density(A_avg)
A_sub_sorted, A_sub_density = get_density(A_sub)
S_avg_sorted, S_avg_density = get_density(S_avg)
S_sub_sorted, S_sub_density = get_density(S_sub)

idx = np.argmax(A_avg_density)
print(A_avg_sorted[idx])

plt.hist(S_avg, bins=20, normed=True, alpha=0.5)
plt.hist(S_sub, bins=20, normed=True, alpha=0.5)
plt.plot(S_avg_sorted, S_avg_density)
plt.plot(S_sub_sorted, S_sub_density)
plt.show()

plt.hist(A_avg, bins=20, normed=True, alpha=0.5)
plt.hist(A_sub, bins=20, normed=True, alpha=0.5)
plt.plot(A_avg_sorted, A_avg_density)
plt.plot(A_sub_sorted, A_sub_density)
plt.show()
