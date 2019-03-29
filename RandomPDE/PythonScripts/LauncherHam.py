import os
import numpy as np

# Integration parameters
hVec = np.zeros(1)
hVec[0] = 0.2
for i in range(1, hVec.size):
    hVec[i] = hVec[i-1] / 2
print(hVec)

M = 200
T = 10000
nT = 1
p = 2.5
det = " "
all = " -all "
program = "Hamiltonian"

for h in hVec:

    outputFileName = "../ReviewPaper/plotStdDev"# + str(h * 1000) #+ "det"
    options = "-h "+ str(h) + " -T " + str(T) + det + " -nT " + str(nT) + \
        " -M " + str(M) + " -p " + str(p) + " -output " + outputFileName + all

    command = "../build/RTSRK-Test/" + program + " " + options

    print(command)
    os.system(command)
