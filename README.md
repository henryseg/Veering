# Veering
Regina-python and sage code for working with transverse taut and veering ideal triangulations. Implemented by Anna Parlak, Saul Schleimer, and Henry Segerman. 

### Download and testing

You will need regina installed inside of python (and/or sage) for the code to run.  Installation instructions can be found here:

https://regina-normal.github.io/  
http://sageregina.unhyperbolic.org/

After cloning or downloading the code, move into the new directory (called 'Veering' if cloned) and check the test suite using

    python test_suite.py
    sage -python test_suite.py

### Usage

There is a list of all isomorphism signatures, up to 16 tetrahedra, in the file Data/veering_census.txt.  As an example of usage, to compute the "big" veering polynomial (defined by Sam Taylor and Michael Landry) start a sage session in the directory and run

    sage: import veering_polynomial
    sage: veering_polynomial.big_polynomial('cPcbbbiht_12')

### Webpage

For references, for information about the census, and for many diagrams, please see:

https://math.okstate.edu/people/segerman/veering.html
