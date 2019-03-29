import os

program = "InvDens"
outFile = "../InvariantDensity/InvDensity"

h = 0.001
T = 100
M = 10
nExp = 5

os.system("cd ../build && make -j32 " + program)

optionsTraj = " -h " + str(h) + " -M " + str(M) + \
   " -T " + str(T) + " -output " + outFile + " -nExp " + str(nExp)

commandTraj = "../build/RTSRK-Test/InvDens" + optionsTraj

os.system(commandTraj)