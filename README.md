# Veering

Regina-python and sage code for working with transverse taut and
veering ideal triangulations; implemented by Anna Parlak, Saul
Schleimer, and Henry Segerman.  The "big" and "small" veering
polynomials are defined by Sam Taylor and Michael Landry.  We thank
Nathan Dunfield for many helpful comments (and for some code).

### Installing regina

Essentially all of the veering code relies on regina; some of it
relies on sage.  Installation instructions for regina into python and
into sage can be found here:

http://regina-normal.github.io/
http://sageregina.unhyperbolic.org/

### Download and testing

After cloning or downloading the veering code, move into the new
directory (called 'Veering' if cloned) and check the test suite using:

    python test_suite.py
    sage -python test_suite.py

### Usage

As an example, start a sage session in the directory and type:

    sage: from file_io import parse_data_file
    sage: veering_isosigs = parse_data_file('Data/veering_census.txt')

The list veering_isosigs now contains all taut isomorphism signatures
for the veering triangulations with at most 16 tetrahedra.  These are
ordered lexicographically.

    sage: sig = veering_isosigs[1]; sig
    'cPcbbbiht_12'

This is the taut isomorphism signature for the only known veering
structure on the figure eight knot complement.  The string before the
underscore is the isomorphism signature for the triangulation; the
string after the underscore records the positions of the edges with
dihedral angle pi in each tetrahedron.

    sage: import taut_polytope
    sage: taut_polytope.is_layered(sig)
    True

This taut structure is layered; we deduce that the figure eight knot
is a fibered knot.

    sage: import veering_polynomial
    sage: veering_polynomial.big_polynomial(sig)
    a^3 - 4*a^2 + 4*a - 1
    sage: veering_polynomial.small_polynomial(sig)
    a^2 - 3*a + 1
    sage: veering_polynomial.small_polynomial(sig, mode = 'alexander')
    a^2 - 3*a + 1

Note that the small veering polynomial divides the big veering
polynomial; this is true in general.  For the figure eight knot the
small polynomial equals the Alexander polynomial; this is not true in
general.

### Webpage

For references, for information about the census, and for many
diagrams, please see:

https://math.okstate.edu/people/segerman/veering.html

### Licence

This work is in the public domain.  See the LICENCE for details.
