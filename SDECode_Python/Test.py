from EM import *
from ParEst import *
import matplotlib.pyplot as mplot

OU = OrnUhl()
EM = EulerMaruyama(OU)
p = np.array([1.0, 0.8])
IC = 0.0
T = 10000.0
N = 100000
Y = EM.solve(IC, N, T, p)

mplot.figure()
mplot.plot(Y)
mplot.show()

h = T / N
parEst = ParEst(Y, OU.gradv, h)
A = parEst.drift()
S = parEst.diffusion()

print([A, S])
