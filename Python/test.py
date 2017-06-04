import numpy as np

sum = 0
for i in np.linspace(1, 100, 100):
    sum = sum + i

v = np.array([2, 5, 10])
w = np.array([1, 1, 1])
dotProd = v.dot(w)
print(dotProd)