The code in this directory is for the Spoof Quasiperfect Factorizations project, which is joint work with Jonathon Hales

See [arxiv]

There are five .py files in this repository, all of which are compatible with Python 3.8.5. We outline them below.

quasiPerfectSearchTools.py is the heart of this respository. It contains an implementation of our algorithm to generate all odd spoof quasiperfect numbers with k factors. 

RUNME.py is a simple wrapper for quasiPerfectSearchTools.py. This code will ask you how many factors you want your spoof quasiperfect factorizations to have, and then generates all odd spoof quasiperfect factorizations with that many factors. It also tells you how long your computations took.

The remaining three files break up larger computations in a way that is hopefully more attainable. Given integers j < k, the python file jofkFactors.py generates all j-factor partial factorizations which plausibly might be extended to a k factor spoof factorization, and saves its output to jofkFactors.npy, where j and k are the given integers. Then, kFactorsFromj.py accepts integers j and k (which should be the same as those provided earlier), as well as a starting index i, and proceeds to check from i up through the end of the document, unless its process is interrupted. It saves any spoofs it discovers to quasiPerfectSolutionsWithkFactors.txt (where k is given), and reports the index of the last factorization it has checked at checkedFromi.txt. If the process is interrupted and then restarted, it proceeds from the last index recorded in checkedFromi.txt, rather than from the given starting index i. kFactorsFromjBackwards.py is identical to kFactorsFromj.py, save that it proceeds backwards instead of forwards from the given index, and it saves its progress in the text document checkedBackFromi.txt. 
