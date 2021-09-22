#import math
#import numpy as np
import time
from quasiPerfectSearchTools import *

try:
    k = int(input("We will generate all nontrivial odd spoof quasiperfect factorizations with k factors. k = ? "))
except:
    k = 0

if k < 1:
    print("That is not a valid input.")
else:
    t = time.time()
    print("Here is a list of all nontrivial odd spoof quasiperfect factorizations with k factors:")
    print(spoofSearch(k))
    print("It took " + str(time.time() - t) + " seconds to calculate them all and report them to you.")
