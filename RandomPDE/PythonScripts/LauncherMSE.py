import os
import numpy as np

# Integration parameters

T = 1
nExp = 5
h = 0.125
nMC = 1000
p = 1.5
nExp = 8
nMCExt = 100

# Compile?
make = True

if make:
    os.system("cd ../build && make -j32")

outputFileName = "outputMSE_ETT"

options = "-h "+ str(h) +  \
      " -T " + str(T) + \
      " -nMC " + str(nMC) + \
      " -output " + outputFileName + \
      " -p " + str(p) + " -nExp " + str(nExp) + \
      " -nMCExt " + str(nMCExt)

program = "MSE"

command = "../build/RTSRK-Test/" + program + " " + options

os.system(command)
print(command)