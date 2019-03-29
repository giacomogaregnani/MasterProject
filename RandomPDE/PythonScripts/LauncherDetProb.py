import os

h = 0.001

M = 20
T = 1000
nT = 15
p = 2.5
det = " "

outputFileName = "../ReviewPaper/testDetProb"

options = "-h "+ str(h) + " -T " + str(T) + det + \
    " -nMC " + str(M) + " -p " + str(p) + " -output " + outputFileName

program = "DetProb"

command = "../build/RTSRK-Test/" + program + " " + options

print(command)
os.system(command)