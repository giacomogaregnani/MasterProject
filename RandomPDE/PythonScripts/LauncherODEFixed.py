import os
import numpy as np

# Integration parameters
hVec = np.zeros(1)

hVec[0] = 0.5
for i in range(1, hVec.size):
    hVec[i] = hVec[i-1] / 10
print(hVec)

hText = '05'

# Observations
noise = .1
noiseText = '1'

# MCMC
nMCMC = 100000
propVar = .5
nObs = 1
prob = " "
prob2 = " -prob2 "
nMC = 10
p = 1.5

# Compile?
make = True
program = "InverseODEFixed"
obsFile = "observations1"

if make:
    os.system("cd ../build && make -j32 " + program)

i = 0

for h in hVec:

    outputFileName = "RESULTSASSYR/MCMC"  + "_h_" + hText + "_noise_" + noiseText + "_rts2"

    options = " -h "+ str(h) + " -noise " + str(noise) + " -nMCMC " + str(nMCMC) + \
          " -propVar " + str(propVar) + \
          " -nObs " + str(nObs) + " -p " + str(p) + \
          prob + " -outputFile " + outputFileName + " -obsFile " + obsFile + " -nMC " + str(nMC) + \
          prob2

    command = "../build/InverseODE/" + program + " " + options

    print(command)
    os.system(command)

    i = i + 1