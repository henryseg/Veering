# veering

Python code (using regina, snappy, and sage) for working with transverse taut
and veering ideal triangulations -- implemented by Anna Parlak, Saul Schleimer,
and Henry Segerman.  We thank Nathan Dunfield for many helpful comments (and
for some code).

### Installation

To install (or uninstall) veering inside Sage via the command line, using pip, type:

    sage -pip install veering

or

    sage -pip uninstall veering

For installation into your system's Python, replace `sage -pip` by `pip3`.
Note that the github repository of veering contains further data and scripts
which are not installed by pip.

Almost all of the veering code relies on regina; some of it relies on
snappy and some relies on SageMath.  Other parts rely on the Python
vector graphics package pyx.  Installation instructions for SageMath,
snappy, and regina can be found at the following webpages:

https://doc.sagemath.org/html/en/installation/ \
https://snappy.math.uic.edu/installing.html \
https://github.com/3-manifolds/regina_wheels

### Testing

After installation start a sage session and run the following:

    sage: import veering
    sage: from veering import test_suite
    sage: test_suite.run_tests()

Each test should take at most a few seconds.

### Usage

As a simple example:

    sage: census = veering.veering_census(); len(census)
    87047

The veering census contains the 87047 taut isomorphism signatures
of the veering triangulations with at most 16 tetrahedra.  These are
ordered lexicographically.

    sage: sig = census[1]; sig
    'cPcbbbiht_12'

This is the taut isomorphism signature for the only veering
structure on the figure eight knot complement.  The string before the
underscore is the isomorphism signature for the triangulation; the
string after the underscore records, for each tetrahedron, which two
edges have dihedral angle pi; the other four edges have dihedral angle 
zero.

    sage: from veering import taut_polytope
    sage: taut_polytope.is_layered(sig)
    True

This taut structure is layered; thus the figure-eight knot is fibered.

    sage: from veering import taut_polynomial
    sage: taut_polynomial.taut_polynomial_via_tree(sig)
    a^2 - 3*a + 1
    sage: taut_polynomial.taut_polynomial_via_tree(sig, mode = 'alexander')
    a^2 - 3*a + 1
    sage: from veering import veering_polynomial
    sage: veering_polynomial.veering_polynomial(sig)
    a^3 - 4*a^2 + 4*a - 1

The taut and veering polynomials are defined by Michael Landry, Yair 
Minsky and Sam Taylor.  Note that the taut polynomial divides the veering 
polynomial; this is true in general.  The taut polynomial of this veering 
triangulation is equal to the Alexander polynomial of the underlying 
manifold; this is not true in general.

    sage: sig = census[257]; sig
    'iLLLQPcbeegefhhhhhhahahha_01110221'
    sage: taut_polytope.cone_in_homology(sig)
    [N(1, -1), N(1, 1)]
    
The cone of homology classes carried by the veering triangulation
`iLLLQPcbeegefhhhhhhahahha_01110221` is spanned by the rays passing 
through (1,-1) and (1,1). Landry, Minsky, and Taylor proved that, 
if nonempty, this cone is equal to a cone on a (not necessarily top-dimensional) 
face of the Thurston norm ball. The chosen basis on H^1 is dual to the 
basis of H_1 used to compute the taut and veering polynomials.

### Webpage

For references, for information about the census, and for many diagrams, 
please see:

https://math.okstate.edu/people/segerman/veering.html

### Citation

When citing the census, please use a version of following (updating the 
version number and the year).
```
@Misc{VeeringCensus,
        author = {Giannopoulos, Andreas and Schleimer, Saul and Segerman, Henry},
        title = {A census of veering structures},
        howpublished = {\url{https://math.okstate.edu/people/segerman/veering.html} {YYYY/MM/DD}}
}
```

When citing the codebase, please use a version of the following (updating the 
version number and the year). 
```
@Misc{Veering,
    author = {Anna Parlak and Saul Schleimer and Henry Segerman},
    title = {veering x.y, code for studying taut and veering ideal triangulations},
    howpublished = {\url{https://github.com/henryseg/Veering}},
    year = {20zz},
}
```

The census and/or code base have been cited in the following works. 

Dynamics of veering triangulations: infinitesimal components of their
flow graphs and applications
- Ian Agol, Chi Cheuk Tsang

Flows, growth rates, and the veering polynomial
- Michael P. Landry, Yair N. Minsky, And Samuel J. Taylor

A polynomial invariant for veering triangulations
- Michael Landry, Yair N. Minsky, And Samuel J. Taylor

Endperiodic maps, splitting sequences, and branched surfaces
- Michael P. Landry And Chi Cheuk Tsang

Geometric triangulations of a family of hyperbolic 3–braids
- Barbara Nimershiem

Computation of the taut, the veering and the Teichmüller polynomials
- Anna Parlak

Arbitrarily large veering triangulations with a vanishing taut polynomial
- Anna Parlak

The taut polynomial and the Alexander polynomial
- Anna Parlak

Obstructing Anosov flows on cusped 3-manifolds
- Misha Schmalian

From veering triangulations to link spaces and back again
- Saul Schleimer and Henry Segerman

Veering branched surfaces, surgeries, and geodesic flows
- Chi Cheuk Tsang

Veering triangulations and transverse foliations
- Jonathan Zung

### Contact

Please contact us with any and all suggestions, questions, and/or corrections.

### Licence

This work is in the public domain. See the LICENCE for details.
