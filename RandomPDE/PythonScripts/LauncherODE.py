import os
import numpy as np

# Integration parameters
hVec = np.zeros(1)

hVec[0] = 0.1
for i in range(1, hVec.size):
    hVec[i] = hVec[0] * pow(2, -i)
print(hVec)

#hText = ['02', '01', '005', '0025']
hText = ['small']

# Observations
noise = .001

# MCMC
nMCMC = 200000
propVar = .01
nObs = 1
prob = " -prob "
noisy = " "
nMC = 10
p = 1.5
isGauss = " "

# Compile?
make = True
program = "InverseODEIC"
obsFile = "observationsShort"

if make:
    os.system("cd ../build && make -j32 " + program)

i = 0

for h in hVec:

    outputFileName = "ODE_IC_PROB_TEST" # + "_h_" + hText[i]
    i = i + 1

    options = " -h "+ str(h) + " -noise " + str(noise) + " -nMCMC " + str(nMCMC) + \
          " -propVar " + str(propVar) + \
          " -nObs " + str(nObs) + " -p " + str(p) + \
          prob + " -nMC " + str(nMC) + \
          " -outputFile " + outputFileName + " -obsFile " + obsFile + \
          noisy + isGauss

    command = "../build/InverseODE/" + program + " " + options

    print(command)
    os.system(command)

    if i == 1:
        nMC = 30
    if i == 2:
        nMC = 5