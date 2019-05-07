import os

# Observations
noise = 0.00001
T = 0.5
nObs = 1

make = True
program = "GenObs"
obsFile = "observationsreview"

if make:
    os.system("cd ../build && make -j32 " + program)

optionsObs = " -h " + str(0.0001) + " -noise " + str(noise) + \
    " -nObs " + str(nObs) + " -T " + str(T) + " -outputFile " + obsFile

commandObs = "../build/InverseODE/GenObs" + optionsObs
os.system(commandObs)
print(commandObs)
