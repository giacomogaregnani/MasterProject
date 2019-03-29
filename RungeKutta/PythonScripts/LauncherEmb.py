import os
import numpy as np

# Integration parameters
hVec = np.zeros(1)
hVec[0] = 0.01

for i in range(1, hVec.size):
    hVec[i] = hVec[0] * pow(2, -i)

print(hVec)

nMCMC = 200000
RAM = " "
propVar = 0.05
obsTime = 2
obsNoise = 1e-4

# Compile?
make = False

if make:
    os.system("cd ../build && make -j32")

i = 0

for h in hVec:

    i = i + 1

    outputFileName = "Results" + str(i)

    options = "-h "+ str(h) + " -nMCMC " + str(nMCMC) + \
          RAM + " -prop " + str(propVar) + \
          " -obsNoise " + str(obsNoise) + \
          " -obsTime " + str(obsTime) + \
          " -outFile " + outputFileName

    program = "OdeInverseEmb"

    command = "../build/Test/" + program + " " + options

    os.system(command)
    print(command)