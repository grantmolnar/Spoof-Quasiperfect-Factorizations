#import math
import numpy as np
import time
from quasiPerfectSearchTools import *

#With this document, we take a numpy file containing all candidate j factor partial factorizations for a k factor quasiperfect spoof, and we compute all 6 factor quasiperfect spoofs

#The number of factors we have
j = int(input("j = ? "))
#The number of factors we are aiming for
k = int(input("k = ? "))
#How often we update how far we've gone
N = 1

#This is the file where our j factor candidates live
fileInput = str(j) + "of" + str(k) + "Factors.npy"
partialFactorsArray = np.load(fileInput)

#This is the number of k-candidate factors our code generated with j factors
length = partialFactorsArray.shape[0]

#We will save any
fileSolutions = "quasiPerfectSolutionsWith" + str(k) + "Factors.txt"

startFromHere = input("Where do you want to start checking? ")

checkedThrough = "checkedBackFrom" + startFromHere + ".txt"

startFromHere = int(startFromHere)

#If you've run code starting from startFromHere before, we'll resume from where you left off
try:
    g = open(checkedThrough, "r")
    startingIndex = int(g.read())
    g.close()
#Otherwise, we'll start from startFromHere
except:
    startingIndex = startFromHere

#We work backwards from our starting index towards 0
for i in range(0, startingIndex + 1):
    spoof = partialFactorsArray[startingIndex-i].tolist()
#    print(spoof)
    #For each candidate partial factorization, we generate all nontrivial k-factor spoof OQFs compatible with that partial factorization
    solutions = spoofFromSeed(k, spoof)
    if len(solutions) > 0:
        #We print our solutions as we find them
        print(solutions)
        #We write solutions, if there are any
        f = open(fileSolutions, "a")
        for factorization in solutions:
            f.write('%s\n' % factorization)
        f.close()
        #Then we make sure we mark that we got that far!
        g = open(checkedThrough, "w")
        g.write(str(startingIndex-i))
        g.close()
    #We inform you every N possibilities we check, and record our progress in checkedThrough
    if i % N == 0:
        print(startingIndex-i)
        g = open(checkedThrough, "w")
        g.write(str(startingIndex-i))
        g.close()
