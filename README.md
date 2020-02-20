# Veering

Regina-python and sage code for working with transverse taut and
veering ideal triangulations.  Implemented by Anna Parlak, Saul
Schleimer, and Henry Segerman.  The "big" and "small" veering
polynomials were first defined by Sam Taylor and Michael Landry.

### Installing regina

All of the veering code relies on regina; some relies on sage.
Installation instructions for regina into python and into sage can be
found here:

http://regina-normal.github.io/
http://sageregina.unhyperbolic.org/

### Download and testing

After cloning or downloading the veering code, move into the new
directory (called 'Veering' if cloned) and check the test suite using

    python test_suite.py
    sage -python test_suite.py

### Usage

There is a list of all isomorphism signatures, up to 16 tetrahedra, in
the file Data/veering_census.txt.  As an example of usage, to compute
the "big" and "small" veering polynomials, start a sage session in the
directory and run

    sage: sig = "cPcbbbiht_12"
    sage: tri, angle = isosig_to_tri_angle(sig)
    sage: veering_polynomial.big_polynomial(tri, angle)
    a^3 - 4*a^2 + 4*a - 1
    sage: veering_polynomial.small_polynomial(tri, angle)
    a^2 - 3*a + 1

### Webpage

For references, for information about the census, and for many
diagrams, please see:

https://math.okstate.edu/people/segerman/veering.html

### Licence

This work is in the public domain.  See the licence for details.
