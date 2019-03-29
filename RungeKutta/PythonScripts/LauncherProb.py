import os
import numpy as np

# Integration parameters
hVec = np.zeros(4)
hVec[0] = 0.1

for i in range(1, hVec.size):
    hVec[i] = hVec[0] * pow(2, -i)

print(hVec)

nMCMC = 100000
RAM = " "
propVar = 0.05
nMC = 800
obsFile = "observations"
noisy = " "

# Compile?
make = False

if make:
    os.system("cd ../build && make -j32")

i = 0

for h in hVec:

    i = i + 1

    outputFileName = "ResultsProb" + str(i)

    options = "-h "+ str(h) + " -nMCMC " + str(nMCMC) + \
          RAM + " -prop " + str(propVar) + \
          " -outFile " + outputFileName + \
          " -nMC " + str(nMC) + noisy

    program = "OdeInverse"

    command = "../build/Test/" + program + " " + options

    os.system(command)
    print(command)