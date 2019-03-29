import os
import numpy as np

# Integration parameters
hVec = np.zeros(1)
hVec[0] = 0.01

for i in range(1, hVec.size):
    hVec[i] = hVec[0] * pow(2, -i)

print(hVec)

p = 2.5
T = 50
nM = 10
nReps = 100

# Compile?
make = False

if make:
    os.system("cd ../build && make -j32")

i = 0
for h in hVec:

    i = i + 1

    for j in range(1, nReps + 1):

        outputFileName = "Measure" + str(i) + "_" + str(j)

        options = "-h "+ str(h) + " -nM " + str(nM) + \
              " -T " + str(T) + " -p " + str(p) + \
              " -output " + outputFileName

        program = "measure"

        command = "../build/Test/" + program + " " + options

        os.system(command)
        print(command)