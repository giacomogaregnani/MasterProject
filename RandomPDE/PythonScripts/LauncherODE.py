import os
import numpy as np

# Integration parameters
hVec = np.zeros(4)

hVec[0] = 0.2
for i in range(1, hVec.size):
    hVec[i] = hVec[0] * pow(2, -i)
print(hVec)

hText = ['02', '01', '005', '0025']

# Observations
noise = .0005

# MCMC
nMCMC = 200000
propVar = .0002
nObs = 1
prob = " -prob "
noisy = " "
nMC = 200
p = 2
isGauss = " -Gauss "
add = " "

# Compile?
make = True
program = "InverseODEIC"
obsFile = "obsreview"

if make:
    os.system("cd ../build && make -j32 " + program)

i = 0

for h in hVec:

    outputFileName = "ODE_IC_PROB_KEPLER_REVIEW" + "_h_" + hText[i]
    i = i + 1

    options = " -h "+ str(h) + " -noise " + str(noise) + " -nMCMC " + str(nMCMC) + \
          " -propVar " + str(propVar) + \
          " -nObs " + str(nObs) + " -p " + str(p) + \
          prob + " -nMC " + str(nMC) + \
          " -outputFile " + outputFileName + " -obsFile " + obsFile + \
          noisy + isGauss + add

    command = "../build/InverseODE/" + program + " " + options

    print(command)
    os.system(command)

    if i == 1:
        nMC = 30
    if i == 2:
        nMC = 5