import math
import copy
import time

# Sieve of erastothenes, returns list of primes up to n which are congruent to 5 or 7 mod 8
def primeSieve(n):
    if n < 2:
        return []
    size = n//3+(n%6==2)
    sieve = [True]*size
    sieve[0] = False
    for i in range(int(n**0.5)//3+1):
        if sieve[i]:
            k=3*i+1|1
            for j in range(k**2//3,size,2*k):
                sieve[j] = False
            for j in range((k*k+4*k-2*k*(i&1))//3,size,2*k):
                sieve[j] = False
    primes = []
    #After this for loop, primes contains all primes from 5 to n
    for j in range(size):
        if sieve[j]:
            primes.append((3*j+1)|1)
    #We now filter for the primes we care about
    return [p for p in primes if p % 8 == 5 or p % 8 == 7]

#primeBound = 1000000

primeBound = int(input("Presieve by primes up to...? "))

#A list of all primes up to primeBound which are congruent to 5 or 7 mod 8
primes5or7mod8 = primeSieve(primeBound)

#Returns True if n is not divisible by any prime in prime5or7mod8, and False otherwise
def notDivisible(n):
    #_iterAllPowers is defined so that n is congruent to 1 or 3 mod 8, so we don't need to worry about n = p == 5 or 7
    #_iterAllPowers also ensures us that n is odd, so the smallest factor of n is at least 3
    bound = n//3
    for p in primes5or7mod8:
        if p <= bound:
            if n % p == 0:
                return False
        else:
            break
    return True

#We think of rational numbers as ordered pairs (x, y) of integers, with (a, b) >= (c, d) if ad - bd > 0

#This is infinity.
#I used to have oo = "Inf", but this is easier on types for numpy
oo = -1

#Given an integer x, returns its sign
sign = lambda x: (x > 0) - (x < 0)

#Given an ordered pair alpha = (alpha0, alpha1) representing a rational number, returns its analogue in reduced terms
def simplify(alpha):
    d = math.gcd(alpha[0], alpha[1])
    return (sign(alpha[1])*alpha[0] // d, sign(alpha[1])*alpha[1] // d)

#If alpha0/alpha1 < beta0/beta1, returns True
#Otherwise, returns False
def lessThan(alpha, beta):
    return alpha[0]*beta[1] < alpha[1]*beta[0]

#If alpha0/alpha1 > beta0/beta1, returns True
#Otherwise, returns False
def greaterThan(alpha, beta):
    return alpha[0]*beta[1] > alpha[1]*beta[0]

#Given a pair alpha = (alpha0, alpha1) representing a rational number, returns alpha0/alpha1 as a float
def toFloat(alpha):
    return alpha[0]/alpha[1]

#Given q and a, returns (q**(a+1)-1)/(q - 1)
def sigma(q, a):
    if q == 1:
        return a + 1
    else:
        return (q**(a+1) - 1)//(q-1)

#Given a spoof factorization [(q1, a1), ..., (qk, ak)], returns tsigma_{-1} of that factorization
def tsigma(spoof):
    #Given a factor (q, a), returns (q - q^-a)/(q - 1)
    def _tsigmaFactor(factor):
        q = factor[0]
        a = factor[1]
        #We don't want trivial factors showing up!
        #assert q != 0 and q != -1
        #If q = 1 and a = oo, we need special syntax, and we need to think extra carefully about this case
        if a == oo:
            return (q, q - 1)
            #return q/(q-1)
        else:
            return (sigma(q,a), q**a)
            #return(q - q**-a)/(q-1)
    numerator = 1
    denominator = 1
    #output = 1
    #We take a product over all factors in our spoof
    for factor in spoof:
        term = _tsigmaFactor(factor)
        #output *= term
        numerator *= term[0]
        denominator *= term[1]
    return simplify((numerator, denominator))
    #return output

#Given a partial factorization partialSpoof, returns a lower bound on tsigma compatible with that partial factorization
#This bound is given subject to the understanding that every term has even power
#This is the function L defined in Theorem 4.4 of [arxiv]
def L(partialSpoof):
    #Given a partial factor, returns the lower extremal associated factor
    #We are assuming that the exponent is even
    def _reduce(partialFactor):
        if partialFactor[0] > 0:
            return (partialFactor[0], partialFactor[1])
        else:
            return (partialFactor[0], partialFactor[2])
    spoof = [_reduce(partialFactor) for partialFactor in partialSpoof]
    return tsigma(spoof)

#Given a partial factorization partialSpoof, returns an upper bound on tsigma compatible with that partial factorization
#This bound is given subject to the understanding that every term has even power
#This is the function U defined in Theorem 4.4 of [arxiv]
def U(partialSpoof):
    #Given a partial factor, returns the upper extremal associated factor
    #We are assuming that the exponent is even
    def _reduce(partialFactor):
        if partialFactor[0] > 0:
            return (partialFactor[0], partialFactor[2])
        else:
            return (partialFactor[0], partialFactor[1])
    spoof = [_reduce(partialFactor) for partialFactor in partialSpoof]
    return tsigma(spoof)

#Given a (partial) spoof factorization, returns the (smallest) number it decomposes.
def eval(partialSpoof):
    output = 1
    for partialFactor in partialSpoof:
        output *= partialFactor[0]**partialFactor[1]
    return output

#Returns what our quasiperfect spoof ought to evaluate to under tsigma_-1
def target(S):
    return (2*S + 1, S)
    #return 2 + 1/S

#Given an odd factor q, yields all powers of r which might occur in an odd spoof quasiperfect factorization
#If lowerBound is given, only returns values greater than or equal to LowerBound
#See Corollary 3.4 of [arxiv]
def _iterAllPowers(q, lowerBound = 0):
    if q % 8 == 1:
        #If our lower bound is less than or equal to 0, we start for 0
        if lowerBound <= 0:
            i = 0
        #If lowerBound is positive and congruent to 0 mod 8, index from lowerBound
        elif lowerBound % 8 == 0 :
            i = lowerBound
            yield i
        #If lowerBound is positive and congruent to a value between 3 and 7 mod 8, index from the first value congruent to 0 mod 8
        elif lowerBound % 8 > 2:
            i = 8*(lowerBound // 8) + 8
            yield i
        #If lowerBound is 1 or 2 mod 8, start from the first value less than lowerBound and congruent to 0 mod 8
        else:
            i = 8*(lowerBound // 8)
        while True:
            yield i + 2
            i += 8
            yield i
    elif q % 8 == 3:
        #If our lower bound is less than or equal to 0, we start for 0
        if lowerBound <= 0:
            i = 0
        #If lowerBound is positive, stat from the first value greater than or equal to lowerBound congruent to 0 mod 4
        elif lowerBound % 4 == 0:
            i = lowerBound
            yield i
        else:
            i = 4*(lowerBound // 4) + 4
            yield i
        while True:
            i += 4
            yield i
    elif q % 8 == 5:
        #If our lower bound is less than or equal to 0, we start for 0
        if lowerBound <= 0:
            i = 0
        #If lowerBound is 7 mod 8, round up to 0 mod 8
        elif lowerBound % 8 == 0:
            i = lowerBound
            yield i
        elif lowerBound % 8 == 7:
            i = 8*(lowerBound // 8) + 8
            yield i
        #If lower bound is between 1 and 6 mod 8, round down to 8
        else:
            i = 8*(lowerBound // 8)
        while True:
            yield i + 6
            i += 8
            yield i
    elif q % 8 == 7:
        #If our lower bound is less than or equal to 0, we start for 0
        if lowerBound <= 0:
            i = 0
        #If lowerBound is positive, stat from the first value greater than or equal to lowerBound congruent to 0 mod 4
        elif lowerBound % 2 == 0:
            i = lowerBound
            yield i
        else:
            i = lowerBound - 1
        while True:
            i += 2
            yield i
    else:
        yield False

#Given a base q, a lowerBound, and an upperBound, yields all (q, a, a) with lowerBound <= a <= upperBound
#such that sigma(q, a) is not divisible by any prime in primes5or7mod8
def _iterFactors(q, lowerBound, upperBound):
    for power in _iterAllPowers(q, lowerBound):
        #We only need to go up to the upperBound
        if upperBound == oo or power < upperBound:
            if notDivisible(sigma(q, power)): #We filter by confirming sigma(q^a) is not divisible by any prime p in primes5or7mod8
                yield (q, power, power)
        else:
            #Once we finish everything up to A (exclusive), we tack on (q, A, oo)
            yield (q, upperBound, oo)
            break

#Given a base r, we compute the the minimal power of r that could occur in our factorization
def minimalPower(r, partialSpoof = []):
    if len(partialSpoof) > 0 and abs(r) == abs(partialSpoof[-1][0]):
        m = partialSpoof[-1][1]
    else:
        m = 2
    for i in _iterAllPowers(r, m):
        return i
        #This break is sort of redundant, but we're really only calling the first entry in iterAllPowers(r)
        break

#Given a list L [[option11, option12, ...], [option21, options22, ...], ...] of lists, returns all lists that take exactly one option from each sublist
#We impose the additional condition that if p = q and p is of lower index than q, then the lower power of p is no greater than the lower power of q
def _product(L, index = 0, output = []):
    if index == 0:
        return _product(L, index = 1, output = [[x] for x in L[0]])
    if index == len(L):
        return output
    else:
        #This function returns True if factor1 and factor2 have unequal bases, or else the lower power of factor1 is less than or equal to the lower power of factor2
        def check(factor1, factor2):
            return factor1[0] != factor2[0] or factor1[1] <= factor2[1]
        #We check here to make sure that if two bases are the same, their next index does not decrease
        return _product(L, index = index + 1, output = [x + [y] for x in output for y in L[index] if check(x[0], y)])

#Given two partial factorizations {(qi, bi, ci)} and {(qi, di, ei)} with the same bases and bi <= di for all i,
#Returns a list of all permissible spoof factorizations {(qi, ai, ai)} with bi <= ai <= di
#such that sigma(qi, ai, ai) is not divisible by any of the primes in prime5or7mod8
def _allCandidates(lowerSpoof, upperSpoof):
    #Our lower and upper spoof factorizations should be the same length
    assert len(lowerSpoof) == len(upperSpoof)
    #For each spoof factor, we compute a list of all possible factors with that base for our spoof satisfying the given bounds
    candidateFactors = []
    for i in range(len(lowerSpoof)):
        #If bi = di, then either bi > B and the factor was left unchanged, or bi = ci and the factor was left unchanged
        if lowerSpoof[i][1] == upperSpoof[i][1]:
            candidateFactors.append([lowerSpoof[i]])
        else:
            #We append (q, a, a), (q, a + 2, a + 2) (or whatever comes next), (q, A - 2, A - 2) (or whatever comes before A), (q, A, oo)
            candidateFactors.append([newFactor for newFactor in _iterFactors(lowerSpoof[i][0], lowerSpoof[i][1], upperSpoof[i][1])])
    return _product(candidateFactors)

#Given a partial spoof and a bound B, replaces all lower bounds on powers which could be infinite with the maximum of B and that factor's lower bound
def _sB(partialSpoof, B):
    #We make it a deep copy as a precaution
    S = copy.deepcopy(partialSpoof)
    for i in range(len(S)):
        if S[i][2] == oo:
            S[i] = (S[i][0], max(S[i][1], B), oo)
    return S

#Given a partialSpoof and a projected number of factors k, returns a disjoint list of partial factorizations that together encode the same spoof factorizations as partialSpoof
#We replace many infinite powers with finite powers, and leave one extra factorization Sa that keeps those infinite powers but raises the lower bound enough to land in case 1 or 2
#and one refinement S with infinite power but L(S) > 2 + 1/S or U(S) < 2
def refineFactorization(partialSpoof, k):
    lowerBound = L(partialSpoof)
    upperBound = U(partialSpoof)
    targetValue = target(eval(partialSpoof))
    #If L(S) <= 2 + 1/S <= U(S), there is hope of a solution!
    if not greaterThan(lowerBound, targetValue) and not greaterThan(targetValue, upperBound):
        #If lowerBound == upperBound and lowerBound <= 2 + 1/S <= upperBound, then we already have a spoof quasiperfect factorization, there is nothing to refine
        if lowerBound == upperBound:
            return [partialSpoof]
        #If lowerBound < upperBound, we will refine our spoof further
        else:
            B = 0
            #We raise B by multiples of 8 until either 2 + 1/S < L(S) or U(S) < 2
            stillInCase = True
            while stillInCase:
                B += 8
                S = _sB(partialSpoof, B)
                lowerBound = L(S)
                upperBound = U(S)
                targetValue = target(eval(S))
                stillInCase = lessThan(lowerBound, upperBound) and not greaterThan(lowerBound, targetValue) and not greaterThan(targetValue, upperBound)
            #We generate all spoofs with indices between the bounds posed by partial spoof and by S
            refinedSpoofs = _allCandidates(partialSpoof, S)
            if len(partialSpoof) < k:
                return refinedSpoofs
            #We know that spoof factorizations with lower bound B on the factors can't work if we can't introduce new factors
            else:
                return refinedSpoofs
    #If all of our factors are included, but lowerBound > targetValue or upperBound < targetValue, there's no hope of fixing things
    elif len(partialSpoof) == k:
        return []
    #Otherwise, just pass the argument through
    else:
        return [partialSpoof]


#Given a (partial) spoof [[q1, b1, c1], ..., [qr, br, cr]], returns the minimum qminus and maximum qplus values for q which are at most and at least -3 and +3 respectively
def _extremalBase(partialSpoof):
    qminus = -3
    qplus = 3
    for factor in partialSpoof:
        qminus = min(qminus, factor[0])
        qplus = max(qplus, factor[0])
    return (qminus, qplus)

#Given a spoof factorization, returns all possible spoof factorizations with one more factor that might work
#In this definition, we presume that any factors of the form 1^b have already been added to our function
#We do not add any trivial factors
def nextFactor(partialSpoof, k):
    #m  is the number of factors in spoof
    m = len(partialSpoof)
    #If our partial spoof already has k factors, we can't add any factors
    if m >= k:
        return []
    lowerBound = L(partialSpoof)
    upperBound = U(partialSpoof)
    spoofValue = eval(partialSpoof)
    #These numbers provide lower bounds on the next factor in our factorization.
    #We presume that +/- 3 are the lowest possible bases, other than 1. If we could increase that bound, then we could increase our lower bound on q
    qminus, qplus = _extremalBase(partialSpoof)
    #Case 1: L(spoof) > 2 + 1/eval(spoof)
    if greaterThan(lowerBound, target(spoofValue)):
        boundOnBase = 1 + math.floor(1/(((2 + 1/spoofValue)/toFloat(lowerBound))**(1/(k-m))-1))
        #Given a candidate factor r^a, verifies that r is not so negative that no spoof perfect number can exist
        def _testU(candidateFactor):
            #The least that our factor of r can reduce
            minimumReduction = L([candidateFactor])
            #The biggest increase possible is if each new base is qplus
            maximumIncreasePerFactor = tsigma([(qplus,oo)])
            numeratorRHS = upperBound[0]*minimumReduction[0]*maximumIncreasePerFactor[0]**(k-m-1)
            denominatorRHS = upperBound[1]*minimumReduction[1]*maximumIncreasePerFactor[1]**(k-m-1)
            return lessThan((2,1),(numeratorRHS, denominatorRHS))
        output = []
        for r in range(qminus,boundOnBase,-2):
            #These powers are not wrong, but are they optimal?
            candidateFactor = (r,minimalPower(r),oo)
            #These tests appear to destroy all nontrivial cases for k = 2 and k = 3. Is this a good thing, or an error?
            if _testU(candidateFactor):
                output.append(partialSpoof + [candidateFactor])
        return output
        #return [partialSpoof + [(r,minimalPower(r),oo)] for r in range(qminus,boundOnBase,-2)]
    #Case 2: U(spoof) <= 2
    elif not greaterThan(upperBound, (2,1)):
        boundOnBase = 1 + math.ceil(1/((2/toFloat(upperBound))**(1/(k-m))-1))
        #Given a candidate exponent r, verifies that r is not so positive that no spoof perfect number can exist
        def _testL(candidateFactor):
            r = candidateFactor[0]
            #The least that our factor of r can reduce
            maximumIncrease = U([candidateFactor])
            #The most extremal decrease possible is if each remaining factor is qminus
            minimumDecreasePerFactor = tsigma([(qminus, oo)])
            #The smallest value Q could take
            boundOnQ = spoofValue*min(r,-qminus)**(2*(k-m))
            boundLHS = target(boundOnQ)
            numeratorRHS = lowerBound[0]*maximumIncrease[0]*minimumDecreasePerFactor[0]**(k-m-1)
            denominatorRHS = lowerBound[1]*maximumIncrease[1]*minimumDecreasePerFactor[1]**(k-m-1)
            #verifies 2 + 1/S > RHS
            return greaterThan(boundLHS,(numeratorRHS, denominatorRHS))
        #This list contains all our candidate partial spoofs after one iterations
        output = []
        for r in range(qplus,boundOnBase,2):
            #We determine, in the weakest possible terms, our next partial factor with base r
            candidateFactor = (r,minimalPower(r),oo)
            #These tests appear to destroy all nontrivial cases for k = 2 and k = 3. Is this a good thing, or an error?
            if _testL(candidateFactor):
                output.append(partialSpoof + [candidateFactor])
        return output
        #return [partialSpoof + [(r,minimalPower(r),oo)] for r in range(qplus,boundOnBase,2)]
    #Case 3a: L(spoof) = U(spoof)
    elif lowerBound == upperBound:
        #Given a candidate exponent r, verifies that r is not so negative that no spoof perfect number can exist
        boundOnBase = 1 + math.floor(1/((50/(98*spoofValue + 49))**(1/(k-m))-1))
        def _testU(candidateFactor):
            r = candidateFactor[0]
            #The least that our factor of r can reduce
            minimumReduction = L([candidateFactor])
            #Recall that we put positive terms before negative terms, so the biggest increase possible is if each new base is abs(r) + 2
            maximumIncreasePerFactor = tsigma([(qplus,oo)])
            numeratorRHS = upperBound[0]*minimumReduction[0]*maximumIncreasePerFactor[0]**(k-m-1)
            denominatorRHS = upperBound[1]*minimumReduction[1]*maximumIncreasePerFactor[1]**(k-m-1)
            #verifies 2 < RHS
            return lessThan((2,1),(numeratorRHS, denominatorRHS))
        output = []
        for r in range(qminus,boundOnBase,-2):
            #We determine, in the weakest possible terms, our next partial factor with base r
            candidateFactor = (r,minimalPower(r),oo)
            if _testU(candidateFactor):
                output.append(partialSpoof + [candidateFactor])
        return output
        #return [partialSpoof + [(r,minimalPower(r),oo)] for r in range(qminus,boundOnBase,-2)]
    #Case 3b: L(spoof) < U(spoof)
    #We do apply refineFactorization before plugging into this function, but it may need to be applied several times before we get what we need
    else:
        return [newSpoof for refinedSpoof in refineFactorization(partialSpoof, k) for newSpoof in nextFactor(refinedSpoof, k)]
    return 0

#Given a collection of candidate spoof quasiperfect factorizations, checks which ones are in fact quasiperfect
def testCandidates(candidateList):
    output = []
    for candidate in candidateList:
        if target(eval(candidate)) == tsigma(candidate):
            output.append(candidate)
    return output

#Returns all candidate partial spoofs which will have k factors, such that all the factors have base 1
def initializeOnes(k):
    candidateOnes = [[]]
    output = []
    #If a is a power of 1 occuring in the rest of our spoof, then a is at least _powerofOneLowerBound(partialSpoof)
    def _powerofOneLowerBound(partialSpoof):
        if len(partialSpoof) == 0:
            return 2
        else:
            #If this is not the first occurence of 1, WLOG the power of 1 shouldn't get smaller as we go up
            return partialSpoof[-1][1]
    #If a is a power of 1 occuring in the rest of our spoof, then a is less than _powerofOneUpperBound(partialSpoof)
    def _powerofOneUpperBound(partialSpoof):
        return math.floor(3/toFloat(tsigma(partialSpoof))*(4/3)**(k-len(partialSpoof)-1))
    #returns True if it is impossible to augment partialSpoof = {(1,b,b)} with more ones, and false otherwise
    def _enoughOnes(partialSpoof):
        #m is the number of factors we already have
        m = len(partialSpoof)
        #a is the smallest power of 1 permissible for the next factor in a partial Spoof
        a = _powerofOneLowerBound(partialSpoof)
        #If sigma(S)*(a+1) > 3*sigma(-3^oo)^(k-m-1), then we can't add more one's
        #If we can prove that every base other than 1^* is larger than 3 is absolute value, we can improve this bound
        return greaterThan(tsigma(partialSpoof),(3*4**(k-m-1), (a+1)*3**(k-m-1)))
    while candidateOnes != []:
        output += candidateOnes
        #We only append those factors [1, b, b] small enough to be possible, with b + 1 not to be divisible by (small) primes congruent to 5 or 7 mod 8
        candidateOnes = [partialSpoof + [(1,b,b)] for partialSpoof in candidateOnes for b in range(_powerofOneLowerBound(partialSpoof), _powerofOneUpperBound(partialSpoof), 2) if notDivisible(b+1)]
    return output

#Given a factor k, and a partial spoof, determines all spoof quasiperfect numbers with k bases that begin with the partial spoof
def spoofFromSeed(k, partialSpoof = []):
    #If we have too many factors, we can't have a k-factor spoof
    if len(partialSpoof) > k:
        return []
    #If we already have all our factors, we just need to pin down what their powers can be
    elif len(partialSpoof) == k:
        return refineFactorization(partialSpoof, k)
    #If we have too few factors, we need to append some more
    else:
        candidateList = [partialSpoof]
        #We are only running k - len(partialSpoof) times, so it is important to make sure every step of our code adds a factor
        for i in range(len(partialSpoof),k):
            #Generate the next factor
            candidateList = [newSpoof for partialSpoof in candidateList for newSpoof in nextFactor(partialSpoof, k)]
            #Refine if necessary.
            candidateList = [refinedSpoof for partialSpoof in candidateList for refinedSpoof in refineFactorization(partialSpoof, k)]
    return testCandidates(candidateList)

#Given the number of factors k, returns all odd spoof quasiperfect factorizations with k factors
def spoofSearch(k):
    candidateList = initializeOnes(k)
    output = []
    for candidate in candidateList:
        output += spoofFromSeed(k, partialSpoof = candidate)
    #Out outputs still need to be checked to see if they are spoof quasiperfect
    return output

#k equals number of bases
#given a spoof and a number of bases, writes that spoof to the relevant document
def writeSpoof(spoof, k):
    assert k > 0
    if k == 1:
        filename = "quasiPerfectWith1Factor.txt"
    else:
        filename = "quasiPerfectWith" + str(k) +"Factors.txt"
    file = open(filename, "a")
    print(spoof)
    file.write(str(spoof))
    file.write("\n")
    file.close()
    return 0

#Given a number of factors k, saves all nontrivial odd spoof quasiperfect numbers with k factors
def spoofReport(k):
    spoofList = spoofSearch(k)
    for spoof in spoofList:
        writeSpoof(spoof, k)
    writeSpoof("Done",k)
    return 0
