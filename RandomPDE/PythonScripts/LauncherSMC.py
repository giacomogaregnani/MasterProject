import os
import numpy as np

# Integration parameters
hVec = np.zeros(5)

hVec[0] = 0.1
for i in range(1, hVec.size):
    hVec[i] = hVec[0] * pow(2, -i)
print(hVec)

# Observations
noise = .05

# MCMC
nParticles = 500000
T = 1
nObs = 10

# Compile?
make = False
program = "ODE_SMC"
obsFile = "observations"
prob = " -prob "

if make:
    os.system("cd ../build && make -j32 " + program)

i = 0
for h in hVec:

    i = i + 1
    outputFileName = "SMC" + str(i)

    options = " -h "+ str(h) + " -noise " + str(noise) + " -nPart " + str(nParticles) + \
          " -T " + str(T) + \
          " -nObs " + str(nObs) + \
          " -outputFile " + outputFileName + " -obsFile " + obsFile + prob\

    command = "../build/InverseODE/" + program + " " + options

    os.system(command)
    print(command)