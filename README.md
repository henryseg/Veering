# Veering

Python code for working with transverse taut and veering ideal triangulations; 
implemented by Anna Parlak, Saul Schleimer, and Henry Segerman. The taut and 
veering polynomials are defined by Michael Landry, Yair Minsky and Sam Taylor. 
We thank Nathan Dunfield for many helpful comments (and for some code).

### Installation


To install Veering inside Sage via the command line type:

    sage -pip install git+https://github.com/henryseg/Veering

If you wish to install in Python, replace `sage -pip` by `pip`
in the above command.

Essentially all of the veering code relies on regina; some of it
relies on snappy and some on SageMath. Other parts rely on the 
Python vector graphics package pyx. Installation instructions for 
SageMath, snappy, and regina can be found at the following webpages:

https://doc.sagemath.org/html/en/installation/ \
https://snappy.math.uic.edu/installing.html \
https://github.com/3-manifolds/regina_wheels

### Testing

For a sanity check do, run in a sage console:

    sage: from veering.test_suite import run_tests
    sage: run_tests()

### Usage

As an example, start a sage sessiond and type:

    sage: from veering.file_io import parse_data_file
    sage: veering_isosigs = parse_data_file('veering_census.txt')

The list `veering_isosigs` now contains all taut isomorphism signatures
for the veering triangulations with at most 16 tetrahedra. These are
ordered lexicographically.

    sage: sig = veering_isosigs[1]; sig
    'cPcbbbiht_12'

This is the taut isomorphism signature for the only known veering
structure on the figure eight knot complement. The string before the
underscore is the isomorphism signature for the triangulation; the
string after the underscore records the positions of the edges with
dihedral angle pi in each tetrahedron.

    sage: from veering import taut_polytope
    sage: taut_polytope.is_layered(sig)
    True

This taut structure is layered; we deduce that the figure eight knot
is a fibered knot.

    sage: from veering import taut_polynomial
    sage: taut_polynomial.taut_polynomial_via_tree(sig)
    a^2 - 3*a + 1
    sage: taut_polynomial.taut_polynomial_via_tree(sig, mode = 'alexander')
    a^2 - 3*a + 1
    sage: from veering import veering_polynomial
    sage: veering_polynomial.veering_polynomial(sig)
    a^3 - 4*a^2 + 4*a - 1
    
Note that the taut polynomial divides the veering polynomial; this is 
true in general. The taut polynomial of this veering triangulation is
equal to the Alexander polynomial of the underlying manifold; this is
not true in general.

    sage: sig = veering_isosigs[257]
    sage: from veering import taut_polytope
    sage: taut_polytope.cone_in_homology(sig)
    [N(1, -1), N(1, 1)]
    
The cone of homology classes carried by the veering triangulation 
`veering_isosigs[257]` is spanned by the rays passing through (1,-1) and
(1,1). Landry, Minsky and Taylor proved that, if nonempty, this cone is
equal to a cone on a (not necessarily top-dimensional) face of the Thurston 
norm ball. The chosen basis on H^1 is dual to the basis of H_1 used to 
compute the taut and veering polynomials.

### Webpage

For references, for information about the census, and for many diagrams, 
please see:

https://math.okstate.edu/people/segerman/veering.html

### Citation

When citing the codebase, please use the following (updating the year). 

```
@Misc{Veering,
    author = {Anna Parlak and Saul Schleimer and Henry Segerman},
    title = {Veering, code for studying taut and veering ideal triangulations},
    howpublished = {\url{https://github.com/henryseg/Veering}},
    year = {20xx},
}
```

### Contact

Please do contact us with any and all suggestions, questions, and/or corrections.

### Licence

This work is in the public domain. See the LICENCE for details.
