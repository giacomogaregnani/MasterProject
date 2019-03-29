import os

# Discretization
N  = 500
h = 1.0 / N
# Prior
gammaExp = .1 # scaling constant
deltaExp = .05 # correlation length
discontinuous = " -discontinuous "
#discontinuous = " "

# Observations
nObs = 50
noise = 1e-4

# Markov Chain Monte Carlo
nMCMC = 100000
nEig = 15 #min(50, round(N-1))

propVar = 3e-3

# number of chains for probabilistic computations
nChains = 10

# Compile?
make = True

if make:
    os.system("cd ../build && make -j32")

outputFile = "infDim"
saveSample = " "
#saveSample = " -saveSample "

options = " -h " + str(h) + \
          " -noise " + str(noise) + \
          " -nMCMC " + str(nMCMC) + \
          " -gammaExp " + str(gammaExp) + \
          " -deltaExp " + str(deltaExp) + \
          " -nEig " + str(nEig) + \
          " -nObs " + str(nObs) + \
          " -outputFile " + outputFile + \
          " -propVar " + str(propVar) + \
          " -nChains " + str(nChains) + \
          discontinuous + saveSample

program = "InversePDEMulti"

command = "../build/InversePDE/" + program + " " + options

os.system(command)
print(command)