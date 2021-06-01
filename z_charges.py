#
# z_charges
#

# Goal - generate z_charges, following [BB]

from sage.matrix.constructor import Matrix
from sage.modules.free_module_element import vector

def tet_vector(i, num_tet):
    """
    Gives the tet equation for the i^th tetrahedron.
    """
    out = []
    for j in range(num_tet): 
        if i == j: out.extend([1]*3)
        else: out.extend([0]*3)
    return out

def solution_vector(num_tet, num_cusps):
    """
    The desired sums of the tet, edge, and holonomy equations.
    """
    out = []
    out.extend([1]*num_tet) # pi
    out.extend([2]*num_tet) # 2 pi
    out.extend([0]*2*num_cusps) # zero
    return vector(out)

def sol_and_kernel(M):
    """
    Given a snappy manifold M, returns a solution to the tet, edge,
    and holonomy equations, as well as a basis for the kernel.
    """
    G = M.gluing_equations()
    T = Matrix(tet_vector(i, M.num_tetrahedra()) for i in range(M.num_tetrahedra()))
    E = T.transpose().augment(G.transpose()).transpose()
    S = solution_vector(M.num_tetrahedra(), M.num_cusps())
    return E.solve_right(S), E.right_kernel().basis()
