# Veering
Regina-python and sage code for working with transverse taut and veering ideal triangulations

After cloning or downloading the code, move into the new directory (Veering if cloned) and check the test suite using

python test_suite.py

or

sage -python test_suite.py

There is a list of all isomorphism signatures in the file Data/veering_census.txt.  As an example of usage, to compute the "big" veering polynomial (as defined by Sam Taylor and Micheal Landry) start a sage session in the directory and run

import veering_polynomial

veering_polynomial.big_polynomial('cPcbbbiht_12')
