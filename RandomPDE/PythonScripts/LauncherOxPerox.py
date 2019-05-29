import os

# Observations
T = 200

make = True
program = "OxPerox"
filename = "resultsOxPeroxReview/OxPeroxAdd"
h = 0.05
nMC = 50
method = 1
p = 1

if make:
    os.system("cd ../build && make -j32 " + program)

options = " -h " + str(h) + " -nMC " + str(nMC) + " -method " + str(method) + \
          " -T " + str(T) + " -output " + filename + " -p " + str(p)

commandObs = "../build/RTSRK-Test/OxPerox" + options
os.system(commandObs)
print(commandObs)