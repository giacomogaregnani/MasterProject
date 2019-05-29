import os
import numpy as np

# Integration parameters
#pVec = np.array([2.5, 3.5, 4.5, 5.5]) # This is used for RK4
pVec = np.array([1.5, 2.5, 3.5]) # This is used for ET

T = 1
nExp = 5
h = 0.125
nMC = 1000

# Compile?
make = True

if make:
    os.system("cd ../build && make -j32")

i = 0
for p in pVec:

    i = i + 1
    outputFileName = "outputMSC_ET" + str(i)

    options = "-h "+ str(h) +  \
          " -T " + str(T) + \
          " -nMC " + str(nMC) + \
          " -output " + outputFileName + \
          " -p " + str(p) + " -nExp " + str(nExp)

    program = "MeanSquareConv"

    command = "../build/RTSRK-Test/" + program + " " + options

    os.system(command)
    print(command)