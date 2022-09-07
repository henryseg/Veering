#
# fundamental_group.py
#

from sage.groups.free_group import FreeGroup
from .transverse_taut import edge_side_face_collections
from .fundamental_domain import spanning_dual_tree

@liberal
def fundamental_group(tri, angle, simplified = True):
    """
    Returns the finite presentation of the fundamental group arising
    from the canonical spanning tree.
    
    Input 

    * "tri" -- a regina triangulation.

    * "angle" -- an veering angle structure

    * "simplified" -- returns the Tietze simpification [Havas, 1994]
      of the fundamental group.

    Examples:

        sage: sig = "cPcbbbdxm_10"  
        sage: fundamental_group(sig, False)
        Finitely presented group < x0, x1, x2, x3 | x1^-1*x0^-1*x1^-1*x3*x2*x3, x2^-1*x1^-1*x2^-1*x0*x3*x0, x0 >

    This is the unsimplified presentation of the fundamental group of
    the figure-eight sibling.  Note that we are using the liberal
    decorator here.

        sage: sig = "cPcbbbiht_12"  # figure-eight
        sage: fundamental_group(sig)
        Finitely presented group < x1, x2 | x1^-1*x2^-2*x1^-1*x2*x1*x2^-1*x1*x2 >
    
    This is the simplified presentation of the fundamental group of
    the figure-eight knot complement.
    """
    
    tree = spanning_dual_tree(tri)
    tree_faces = tree[0]
    F = FreeGroup(tri.countTriangles())
    # We add all of the faces as generators.  This is useful, for
    # example, when defining homomorphisms from \pi_1(M).
    gens = F.gens()

    # can this be done with just the routines in fundamental_domain?
    esfc = edge_side_face_collections(tri, angle)

    relators = []
    for edge in esfc:  # relations for the edges
        rel = F.one()
        left, right = edge
        left.reverse()
        for term in left:
            rel = rel * gens[term[0]].inverse()
        for term in right:
            rel = rel * gens[term[0]]
        relators.append(rel)

    gens = F.gens()
    for f in tree_faces:  # killing the tree faces
        relators.append(gens[f])
    G = F/relators
    if simplified:
        G = G.simplified()
    return G
