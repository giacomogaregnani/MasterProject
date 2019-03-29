import os
import numpy as np

# Discretization
NVec = np.zeros(4)
NVec[0] = 20
for i in range(1, NVec.size):
    NVec[i] = NVec[i-1] * 2

# Observations
nObs = 4
noise = 1e-4

# Markov Chain Monte Carlo
nMCMC = 1000
propVar = 3e-2

# number of chains for probabilistic computations
nChains = 200

# output samples
saveSample = " -saveSample "

# Compile?
make = True

if make:
    os.system("cd ../build && make -j32")

for i in range(0, NVec.size):

    N = NVec[i]
    h = 1.0 / N

    if i == 0:
        obsGiven = " "
    else:
        obsGiven = " -obsGiven "

    outputFile = "finDim" + str(int(N))

    options = " -h " + str(h) + \
              " -noise " + str(noise) + \
              " -nMCMC " + str(nMCMC) + \
              " -nObs " + str(nObs) + \
              " -outputFile " + outputFile + \
              " -propVar " + str(propVar) + \
              " -nChains " + str(nChains) + \
              saveSample + \
              obsGiven

    program = "InversePDEFinite"

    command = "../build/InversePDE/" + program + " " + options

    os.system(command)
    print(command)