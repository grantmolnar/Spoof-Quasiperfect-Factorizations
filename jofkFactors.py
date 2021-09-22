#import math
import numpy as np
import time
from quasiPerfectSearchTools import *

#With this file, we generate all k-candidate j-factor partial factorizations suggested by our algorithm

#The number of factors we have
j = int(input("j = ? "))
#The number of factors we are aiming for
k = int(input("k = ? "))

candidateList = [[]]
print(candidateList)

#This code emulates the spoofSearch and spoofFromSeed functions, but stops early
for i in range(1,j+1):
    candidateList = [newSpoof for partialSpoof in candidateList for newSpoof in nextFactor(partialSpoof, k)]
    print(len(candidateList))
    #Refining our factorizations to land in case 1 or 2 of our algorithm can take some time
    #If you set this condition to be True, we would save the partial factorizations we produce prior to that refining process
    if False:
        filename = str(i) + "of" + str(k) + "FactorsUnrefined.npy"
    #    filename = str(i) + "of" + str(k) + "FactorsUnrefined.txt"
        f = open(filename, "wb")
    #    f = open(filename, "w")
        np.save(f, np.array(candidateList))
    #    for factorization in candidateList:
    #        f.write('%s\n' % factorization)
        f.close()
    candidateList = [refinedSpoof for partialSpoof in candidateList for refinedSpoof in refineFactorization(partialSpoof, k)]
    #We report how many candidate spoofs there are after each iteration through a new host of factors
    print(len(candidateList))
    #We could replace the i == j condition with True to save not only our jofkFactors.npy record, but also iofkFactors.npy for all i <= j
    if i == j:
        filename = str(i) + "of" + str(k) + "Factors.npy"
    #    filename = str(i) + "of" + str(k) + "Factors.txt"
        f = open(filename, "wb")
    #    f = open(filename, "w")
        np.save(f, np.array(candidateList))
    #    for factorization in candidateList:
    #        f.write('%s\n' % factorization)
        f.close()
