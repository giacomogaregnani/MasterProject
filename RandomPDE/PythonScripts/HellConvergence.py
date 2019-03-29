import os
import numpy as np

# Integration parameters
nMCVec = np.zeros(1)

nMCVec[0] = 1
for i in range(1, nMCVec.size):
    nMCVec[i] = nMCVec[0] * pow(2, i)
print(nMCVec)

hVec = np.zeros(7)

hVec[0] = 0.1
for i in range(1, hVec.size):
    hVec[i] = hVec[0] * pow(2, -i)
print(hVec)

# Observations
noise = 0.05

# MCMC
nMCMC = 1e6
T = 1
propVar = 0.1
nObs = 10
prob = " -prob "
noisy = " -noisy "

# Compile?
make = True
program = "InverseODE"
obsFile = " observations "

if make:
    os.system("cd ../build && make -j32 " + program)

i = 0

for h in hVec:

    i = i + 1

    # Generate reference distribution (exact algorithm)
    outputFileName = "MCMC_ODE_HELL_REF_h" + str(i)
    options = " -h "+ str(h) + " -noise " + str(noise) + " -nMCMC " + str(nMCMC) + \
              " -propVar " + str(propVar) + " -T " + str(T) + \
              " -nObs " + str(nObs) + \
              prob + " -nMC " + str(1) + \
              " -outputFile " + outputFileName + " -obsFile " + obsFile

    command = "../build/InverseODE/" + program + " " + options
    print(command)
    os.system(command)

    k = 0

    for nMC in nMCVec:

        k = k + 1

        for j in range(0, 1):

            outputFileName = "MCMC_ODE_HELL_h_" + str(i) + "_nMC_" + str(k) + "_" + str(j)

            options = " -h " + str(h) + " -noise " + str(noise) + " -nMCMC " + str(nMCMC) + \
                      " -propVar " + str(propVar) + " -T " + str(T) + \
                      " -nObs " + str(nObs) + \
                      prob + " -nMC " + str(nMC) + \
                      " -outputFile " + outputFileName + " -obsFile " + obsFile + " -noisy "

            command = "../build/InverseODE/" + program + " " + options

            print(command)
            os.system(command)