#
# taut_polytope.py
#

# Analyze manifolds and their taut/veering triangulations.

# Solvers - GLPK, CPLEX, CVXOPT, Gurobi

import regina
import snappy

from sage.numerical.mip import MIPSolverException, MixedIntegerLinearProgram
from sage.geometry.polyhedron.constructor import Polyhedron
from sage.matrix.constructor import Matrix
from sage.modules.free_module_integer import IntegerLattice
from sage.modules.free_module_element import vector
from sage.geometry.cone import Cone
from sage.rings.integer_ring import ZZ

from .taut import liberal
from .transverse_taut import is_transverse_taut
from .taut_homology import edge_equation_matrix_taut, elem_vector, faces_in_smith, rank_of_quotient
from .fundamental_domain import non_tree_face_cycles


# Examining edge/face matrices


# Sage cannot treat the variables of a MILP as a vector, so we do it
# ourselves.
def dot_prod(u, v):
    try:
        dim = len(u)
    except:
        dim = len(v)
    return sum( u[i] * v[i] for i in range(dim) )


def extract_solution(q, v):
    dim = q.number_of_variables()
    return [q.get_values(v[i]) for i in range(dim)]


# unused?
def get_polytope(N):
    """
    Compute the (normalised) polytope cut out by N.
    """
    num_faces = N.dimensions()[1]
    # q = MixedIntegerLinearProgram( maximization = False, solver = "Coin" ) # Why!!!
    q = MixedIntegerLinearProgram( maximization = False, solver = "GLPK" ) # Why!!!
    w = q.new_variable(real = True, nonnegative = True)
    for v in N.rows():
        q.add_constraint( dot_prod(v, w) == 0 )
    q.add_constraint( sum( w[i] for i in range(num_faces) ) == 1 )
    return q.polyhedron()


# never use this
def farkas_solution(N):
    """
    Look for an edge vector u so that all entries of u*N are positive,
    minimizing the sum of the entries of u*N.  If one exists, returns
    (True, u).
    """
    q = MixedIntegerLinearProgram( maximization = False, solver = "GLPK" )
    u = q.new_variable( real = True, nonnegative = False )
    for v in N.columns():
        q.add_constraint( dot_prod(u, v), min = 1 )
    q.set_objective( sum( dot_prod(u, v) for v in N.columns() ) )
    try:
        q.solve()
        U = extract_solution(q, u)
        out = (True, U)
    except MIPSolverException:
        out = (False, None)
    return out


def non_trivial_solution(N, real_bool = True, int_bool = False,
                         solver = "GLPK", upper_bound = None):
    """
    Look for a non-negative, non-zero face vector w with N*w = 0.  If
    one exists, returns (True, w).
    """
    num_faces = N.dimensions()[1]
    # q = MixedIntegerLinearProgram( maximization = False, solver = "Gurobi" )
    # q = MixedIntegerLinearProgram( maximization = False, solver = "CVXOPT" )
    # q = MixedIntegerLinearProgram( maximization = False, solver = "CPLEX" )
    q = MixedIntegerLinearProgram( maximization = False, solver = solver )
    w = q.new_variable ( real = real_bool, integer = int_bool, nonnegative = True )
    # int has topological meaning here - hilbert versus vertex!
    for v in N.rows():
        q.add_constraint( dot_prod(v, w) == 0 )
    S = sum( w[i] for i in range(num_faces) )
    q.add_constraint( S >= 1 )
    if upper_bound != None:
        q.add_constraint( upper_bound >= S )
    # q.set_objective( S )  
    q.set_objective( None )  # any non-zero solution will do - and this is _much_ faster.
    try:
        q.solve()
        W = extract_solution(q, w)
        out = (True, W)
    except MIPSolverException:
        out = (False, None)
    return out


@liberal
def get_non_triv_sol(tri, angle):
    N = edge_equation_matrix_taut(tri, angle)
    N = Matrix(N)
    non_triv, sol = non_trivial_solution(N, real_bool = False, int_bool = True)
    twiddles = is_transverse_taut(tri, angle, return_type = "face_coorientations")
    sol = [int(a * b) for a, b in zip(sol, twiddles)]
    return sol


def fully_carried_solution(N):
    """
    Look for a face vector w with N*w = 0 and with all entries
    positive, minimizing the sum of its entries.  If one exists,
    returns (True, w).
    """
    num_faces = N.dimensions()[1]
    # q = MixedIntegerLinearProgram( maximization = False, solver = "Gurobi" ) # Grrr.
    # q = MixedIntegerLinearProgram( maximization = False, solver = "Coin" ) # Why!!!
    # q = MixedIntegerLinearProgram( maximization = False, solver = "GLPK" )
    q = MixedIntegerLinearProgram( maximization = False, solver = "PPL" )
    w = q.new_variable(real = True, nonnegative = True)
    # w = q.new_variable(integer = True, nonnegative = True)
    for v in N.rows():
        q.add_constraint( dot_prod(v, w) == 0 )
    for i in range(num_faces):
        q.add_constraint( w[i] >= 1 )
    # S = sum( w[i] for i in range(num_faces) ) # Euler char
    # S = q.sum( w[i] for i in range(num_faces) ) # Euler char
    q.set_objective( None )  # any positive solution will do
    try:
        q.solve()
        W = extract_solution(q, w)
        out = (True, W)
    except MIPSolverException:
        out = (False, None)
    return out


@liberal
def is_layered(tri, angle):
    """
    Given a tri, angle decide if it is layered.
    """
    N = edge_equation_matrix_taut(tri, angle)
    N = Matrix(N)
    layered, _ = fully_carried_solution(N)
    return layered

@liberal
def min_carried_neg_euler(tri, angle):
    """
    Given a tri, angle, find the minimal negative euler characteristic
    among carried surfaces.  Returns zero if there is none such.
    """
    N = edge_equation_matrix_taut(tri, angle)
    N = Matrix(N)
    # next line is broken... non_triv need not give the min... 
    exists, sol = non_trivial_solution(N, real_bool = False, int_bool = True)
#    non_trivial_solution(N, real_bool = True, int_bool = False,
#                         solver = "GLPK", upper_bound = None):
    if sol == None:
        out = 0.0
    else:
        out = sum(sol)/2 # sum counts triangles, so divide
    return out


@liberal
def carries_torus_or_sphere(tri, angle):
    """
    Given a tri, angle decides if there is a once-punctured torus or
    three-times punctured sphere among the carried surfaces.
    """
    N = edge_equation_matrix_taut(tri, angle)
    # We must decide if there is a pair of columns which are negatives
    # of each other.  So we use the unique+sort trick.
    N = Matrix(N).transpose() # easier to work with rows
    N = set(tuple(row) for row in N) # make rows unique
    N = Matrix(N)  # make rows vectors
    P = []
    for row in N:
        P.append(row)
        P.append(-row)
    P.sort()
    for i in range(len(P) - 1):
        if P[i] == P[i+1]:
            return True
    return False


@liberal
def is_torus_bundle(tri, angle):
    """
    Given a tri, angle, decides if it is layered by once-punctured
    tori.
    """
    return (tri.homology().rank() == 1 and
            is_layered(tri, angle) and
            carries_torus_or_sphere(tri, angle))


# TODO - check this against the homological dim of the taut
# cone... and then delete.  :)
@liberal
def LMN_tri_angle(tri, angle):
    """
    Given a tri, angle, decide LMN
    """
    N = edge_equation_matrix_taut(tri, angle)
    N = Matrix(N)
    farkas, farkas_sol = farkas_solution(N)
    non_triv, non_triv_sol = non_trivial_solution(N)
    full, full_sol = fully_carried_solution(N)
    if full:
        out = "L"  # depth zero
    elif non_triv:
        out = "M"  # who knows!
    else:
        assert farkas
        out = "N"  # depth infinity 
    return out


@liberal
def analyze_deeply(tri, angle):
    N = edge_equation_matrix_taut(tri, angle)
    N = Matrix(N)
    M = snappy.Manifold(tri.snappystring())

    # look at it
    alex = alex_is_monic(M)
    hyper = hyper_is_monic(M)
    non_triv, non_triv_sol = non_trivial_solution(N)
    full, full_sol = fully_carried_solution(N)
    try:
        assert non_triv or not full # full => non_triv
        assert alex or not full # full => fibered => alex is monic
        assert hyper or not full # full => fibered => hyper is monic
        assert alex or not hyper # hyper is monic => alex is monic
    except AssertError:
        print("contradiction in maths")
        raise
    if full:
        pass
    elif non_triv:
        print("non-triv sol (but not full)")
        print((alex, hyper))
    elif not alex or not hyper:
        print("no sol")
        print((alex, hyper))
    return None


# Code to compute the dimension of the taut cone as projected into
# homology.


# (co)homology


def matrix_transpose(M):
    return list(map(lambda *row: list(row), *M))


@liberal
def zeroth_coboundary(triangulation):
    """
    Given a regina triangulation, returns the zeroth coboundary matrix
    for the dual cell structure
    """
    # for every primal tetrahedron, find its boundary.
    matrix = []
    for tet in triangulation.tetrahedra():
        row = [0] * triangulation.countFaces(2)
        for i in range(4):
            tri_index = tet.triangle(i).index()
            perm = tet.faceMapping(2, i)
            row[tri_index] += perm.sign()
        matrix.append(row)
    return matrix_transpose(matrix)


@liberal
def first_coboundary(triangulation):
    """
    Given a regina triangulation, returns the first coboundary matrix
    for the dual cell structure
    """
    # for each primal face, find the primal edges incident to it.  Put
    # +1/-1 as the orientation around the triangle agrees or disagrees
    # with orientation on this edge
    matrix = []
    for tri in triangulation.triangles():
        row = [0] * triangulation.countEdges() # == countFaces(1)
        for i in range(3):
            edge_index = tri.edge(i).index()
            perm = tri.edgeMapping(i) # == tri.faceMapping(1, i)
            row[edge_index] += perm.sign()
        matrix.append(row)
    return matrix_transpose(matrix)


# extracting the taut cone


def taut_rays_matrix(N):
    # get the extreme rays of the taut cone - note that the returned
    # vectors are non-negative because they "point up"
    dim = N.dimensions()[1]
    elem_ieqs = [[0] + list(elem_vector(i, dim)) for i in range(dim)]
    N_rows = [v.list() for v in N.rows()]
    N_eqns = [[0] + v for v in N_rows]
    P = Polyhedron(ieqs = elem_ieqs, eqns = N_eqns)

    rays = [ray.vector() for ray in P.rays()]
    for ray in rays:
        assert all(a.is_integer() for a in ray)
    # all of the entries are integers, represented in QQ, so we clean them.
    return [vector(ZZ(a) for a in ray) for ray in rays]


@liberal
def taut_rays(tri, angle):
    # get the extreme rays of the taut cone - note that the returned
    # vectors are non-negative because they "point up"
    N = edge_equation_matrix_taut(tri, angle)
    N = Matrix(N)
    return taut_rays_matrix(N)


# the function

# Factor this - produce the image in H_1 with basis the non-tree edges
# and then take dimension.  This is better, because we get a Thurston face
# in the "correct" basis - that is the basis where we computed the
# taut polynomial.  This allows us to correctly compare Newton polytopes.


@liberal
def taut_cone_homological_dim(tri, angle):
    # find the dimension of the projection of the taut cone into
    # homology

    # boundaries of tets
    bdys = zeroth_coboundary(tri)
    bdys = matrix_transpose(bdys)
    rays = taut_rays(tri, angle)
    # but these are all 'upwards', so we need to fix the
    # co-orientations
    coorient = is_transverse_taut(tri, angle, return_type = "face_coorientations")
    rays = [[int(a * b) for a, b in zip(coorient, ray)] for ray in rays]

    # now work in the space of two-chains
    Rays = IntegerLattice(rays + bdys)
    Bdys = Matrix(bdys)
    Cobs = Bdys.transpose()
    Anns = Cobs.kernel()
    return Rays.intersection(Anns).dimension()


@liberal
def projection_to_homology(tri,angle):
    non_tree_as_cycles = non_tree_face_cycles(tri, angle)
    Q = Matrix(non_tree_as_cycles)
    S, U, V = faces_in_smith(tri,angle,[])
    rank, dimU, dimV = rank_of_quotient(S)
    U = Matrix(U)
    P = U.transpose().inverse()
    P = P.delete_rows(range(0, dimU-rank))
    A = P*Q
    return A


@liberal
def cone_in_homology(tri, angle):
    rays = taut_rays(tri, angle)
    if len(rays) == 0:
        return []
    else:
        A = projection_to_homology(tri,angle)
        projectedRays = [A*v for v in rays]
        C = Cone(projectedRays)
        return [ray for ray in C.rays()]


@liberal
def analyse_many_angles(tri):
    angles = regina.AngleStructures.enumerate(tri, True)
    angles = [regina_taut_struct_to_ints(angles.structure(i)) for i in range(angles.size())]
    dict_of_dimensions = {}
    for angle in angles:
        angle_str = "".join(str(a) for a in angle)
        if is_layered(tri, angle):
            if "layered" not in dict_of_dimensions:
                dict_of_dimensions["layered"] = set(angle_str)
            else:
                dict_of_dimensions["layered"].add(angle_str)
        else:
            dim = taut_cone_homological_dim(tri, angle)
            if dim not in dict_of_dimensions:
                dict_of_dimensions[dim] = set(angle_str)
            else:
                dict_of_dimensions[dim].add(angle_str)


# coumputing the depth


def cut_down_matrix(N, u):
    """
    Given N, an edge/face incidence matrix, and u, the sum of the
    extremal rays, return the matrix N' where we (a) delete all face
    columns (for all non-zero entries in u) and (b) delete all edge
    rows (for all non-zero entries in the to-be-deleted face columns).

    The matrix N' is the edge/face incidence matrix for the taut
    triangulation after cutting along the surface given by u (and note
    that we delete product regions).  If N' is empty, then there are
    two possibilities: either there are no internal faces left, or
    there are.  If there are none, then we are done.  If there are
    some, then we have to cut one more time, but there is no matrix
    left - so we return the "error" as a second argument.
    """
    num_edges = N.dimensions()[0]
    num_faces = N.dimensions()[1]
    cols_to_kill = [ i for i in range(num_faces) if u[i] > 0 ]
    rows_to_kill = [ j for j in range(num_edges) if any([  N[j][i] != 0 for i in cols_to_kill  ]) ]

    if len(cols_to_kill) == num_faces: # layered so
        return None, 0

    if len(rows_to_kill) == num_edges: # not layered, but no edges left after cutting so
        return None, 1

    M = [N[j] for j in range(num_edges) if j not in rows_to_kill]
    M = Matrix(M)
    M = M.transpose()
    Np = [M[i] for i in range(num_faces) if i not in cols_to_kill]
    Np = Matrix(Np)
    Np = Np.transpose()
    return Np, 0


@liberal
def depth(tri, angle):
    """
    Given a transverse taut triangulation compute the depth of the
    horizontal branched surface B.  If the depth is finite we return
    (True, depth).  If the depth is not defined then we return (False,
    cuts).  Here cuts is the number of times we must cut before
    getting to a sutured manifold equipped with a veering
    triangulation (possibly with boundary) that carries nothing.
    """
    N = edge_equation_matrix_taut(tri, angle)
    N = Matrix(N)
    cuts = 0
    while True:
        rays = taut_rays_matrix(N)
        if len(rays) == 0:
            return (False, cuts)
        u = sum(rays)
        N, e = cut_down_matrix(N, u)
        if N == None:
            return (True, cuts + e)
        cuts = cuts + 1
