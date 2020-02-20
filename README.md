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

As an example, start a sage session in the directory and type:

    sage: from file_io import parse_data_file
    sage: veering_isosigs = parse_data_file('Data/veering_census.txt')

The list veering_isosigs now contains taut isomorphism signatures for
all veering triangulations, up to 16 tetrahedra.  Our code is mainly
intended for batch processing.  However, examining individual
manifolds is also possible:

    sage: sig = veering_isosigs[1]; sig
    'cPcbbbiht_12'

This is the taut isomorphism signature for the only known veering
structure on the figure eight knot complement.  We convert this to a
regina triangulation and a taut angle structure as follows:

    sage: from taut import isosig_to_tri_angle
    sage: tri, angle = isosig_to_tri_angle(sig)

We can now compute various properties and invariants:

    sage: import taut_polytope
    sage: taut_polytope.is_layered(tri, angle)
    True

So this taut triangulation is layered; thus the figure eight knot
complement is fibred.  To compute the small and large polynomials
type:

    sage: import veering_polynomial
    sage: veering_polynomial.big_polynomial(tri, angle)
    a^3 - 4*a^2 + 4*a - 1
    sage: veering_polynomial.small_polynomial(tri, angle)
    a^2 - 3*a + 1

Note that the small polynomial divides the large.

### Webpage

For references, for information about the census, and for many
diagrams, please see:

https://math.okstate.edu/people/segerman/veering.html

### Licence

This work is in the public domain.  See the licence for details.
