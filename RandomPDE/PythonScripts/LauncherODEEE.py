import os
import numpy as np

# Integration parameters
hVec = np.zeros(1)

hVec[0] = 0.1
for i in range(1, hVec.size):
    hVec[i] = hVec[0] * pow(2, -i)
print(hVec)

hText = ['02', '01', '005', '0025']

# Observations
noise = .001

# MCMC
nMCMC = 200000
propVar = .01
nObs = 1
isGauss = " "

# Compile?
make = True
program = "InverseODEICEE"
obsFile = "observationsShort"
nErr = 200
L = 10

if make:
    os.system("cd ../build && make -j10 " + program)

i = 0

for h in hVec:

    outputFileName = "ODE_IC_ERREST_TEST" # + "_h_" + hText[i]
    i = i + 1

    options = " -h "+ str(h) + " -noise " + str(noise) + " -nMCMC " + str(nMCMC) + \
          " -propVar " + str(propVar) + \
          " -nObs " + str(nObs) + \
          " -nErr " + str(nErr) + \
          " -outputFile " + outputFileName + " -obsFile " + obsFile + \
          isGauss + " -L " + str(L)

    command = "../build/InverseODE/" + program + " " + options

    print(command)
    os.system(command)