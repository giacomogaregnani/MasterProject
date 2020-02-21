import numpy as np
import matplotlib.pyplot as plt

v = np.empty(4)
print v

A = np.zeros((2, 2))
print A

v = np.arange(10) + 1
print v

I = np.eye(2)
print I

a = np.zeros(10)
for i in range(1, 10):
    v = i * np.ones(10)
    iCol = float(i) / 11
    # plt.plot(v, color=(iCol, iCol, iCol))
    plt.plot(v, color='b', alpha=1-iCol)
plt.show()