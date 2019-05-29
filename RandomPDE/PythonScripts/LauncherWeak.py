import os
import numpy as np

# Integration parameters
#pVec = np.array([1.5, 2, 2.5, 3]) # This is used for RK4
pVec = np.array([1, 1.5, 2]) # This is used for ET

T = 1
nExp = 5
h = 0.125
nMC = 1e6

# Compile?
make = True

if make:
    os.system("cd ../build && make -j32")

i = 0
for p in pVec:

    i = i + 1
    outputFileName = "outputWEAK_ET" + str(i)

    options = "-h "+ str(h) +  \
          " -T " + str(T) + \
          " -nMC " + str(nMC) + \
          " -output " + outputFileName + \
          " -p " + str(p) + " -nExp " + str(nExp)

    program = "WeakConv"

    command = "../build/RTSRK-Test/" + program + " " + options

    os.system(command)
    print(command)