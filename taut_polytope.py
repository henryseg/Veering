#
# classify_taut_ideal_triangulations
#

#
# Goal - analyze manifolds and their taut/veering triangulations.
#

# to do - use Gurobi

# for a few example usages see Veering_census/census_reading.py

import regina
import snappy

from sage.numerical.mip import MIPSolverException, MixedIntegerLinearProgram
from sage.geometry.polyhedron.constructor import Polyhedron
from sage.matrix.constructor import Matrix
from sage.modules.free_module_integer import IntegerLattice

from taut import liberal
from transverse_taut import is_transverse_taut
from taut_homology import edge_equation_matrix_taut

# Polynomials that come with snappy (in sage)

def alex_is_monic(M):
    p = M.alexander_polynomial()
    return p.is_monic()

def hyper_is_monic(M): # add a way to dial up the precision
    p = M.hyperbolic_torsion()
    lead = p.coefficients(sparse=False)[-1]
    return abs(1 - lead) < 0.000001  # worry about lead = -1

### Examining edge/face matrices

# Sage cannot treat the variables of a MILP as a vector, so we have to
# do things ourselves.  Boo.
def dot_prod(u, v):
    try:
        dim = len(u)
    except:
        dim = len(v)
    return sum(u[i]*v[i] for i in range(dim))


def extract_solution(q, v):
    dim = q.number_of_variables()
    return [q.get_values(v[i]) for i in range(dim)]


def get_polytope(N):
    """
    Compute the (normalised) polytope cut out by N.
    """
    num_faces = N.dimensions()[1]
    # q = MixedIntegerLinearProgram( maximization = False, solver = 'Coin' ) ### Why!!!
    q = MixedIntegerLinearProgram( maximization = False, solver = 'GLPK' ) ### Why!!!
    w = q.new_variable(real = True, nonnegative = True)
    for v in N.rows():
        q.add_constraint( dot_prod(v, w) == 0 )
    q.add_constraint( sum( w[i] for i in range(num_faces) ) == 1 )
    return q.polyhedron()


def farkas_solution(N):  # never use this
    '''Look for an edge vector u so that all entries of u*N are positive,
    minimizing the sum of the entries of u*N.  If one exists, returns
    (True, u).'''
    q = MixedIntegerLinearProgram( maximization = False, solver = 'GLPK' )
    u = q.new_variable( real = True, nonnegative = False )
    for v in N.columns():
        q.add_constraint( dot_prod(u, v), min = 1 )
    q.set_objective( sum( dot_prod(u, v) for v in N.columns() ) )
    try:
        q.solve()
        U = extract_solution(q,u)
        out = (True, U)
    except MIPSolverException:
        out = (False, None)
    return out


def non_trivial_solution(N, real_bool = True, int_bool = False):
    """
    Look for a face vector w with N*w = 0, minimizing the sum of its
    entries.  If one exists, returns (True, w).
    """
    num_faces = N.dimensions()[1]
    q = MixedIntegerLinearProgram( maximization = False, solver = 'GLPK' )
    w = q.new_variable ( real = real_bool, integer = int_bool, nonnegative = True )
    # actually, int has topological meaning here - hilbert versus vertex!
    # w = q.new_variable ( integer = True, nonnegative = True )
    for v in N.rows():
        q.add_constraint( dot_prod(v, w) == 0 )
    S = sum( w[i] for i in range(num_faces) )
    q.add_constraint( S >= 1 )
    q.set_objective( S )
    try:
        q.solve()
        W = extract_solution(q,w)
        out = (True, W)
    except MIPSolverException:
        out = (False, None)
    return out
        

@liberal
def get_non_triv_sol(tri, angle):
    N = edge_equation_matrix_taut(tri, angle)
    N = Matrix(N)
    non_triv, sol = non_trivial_solution(N, real_bool = False, int_bool = True)
    twiddles = is_transverse_taut(tri, angle, return_type='face_coorientations')
    sol = [int(a*b) for a, b in zip(sol, twiddles)]
    return sol


def vertex_solutions(N, real_bool = True, int_bool = False):
    """
    Look for a face vector w with N*w = 0, minimizing the sum of its
    entries.  If one exists, returns (True, w).
    """
    num_faces = N.dimensions()[1]
    q = MixedIntegerLinearProgram( maximization = False, solver = 'GLPK' )
    w = q.new_variable ( real = real_bool, integer = int_bool, nonnegative = True )
    # actually, int has topological meaning here - hilbert versus vertex!
    # w = q.new_variable ( integer = True, nonnegative = True )
    for v in N.rows():
        q.add_constraint( dot_prod(v, w) == 0 )
    S = sum( w[i] for i in range(num_faces) )
    q.add_constraint( S >= 1 )
    q.set_objective( S )
    try:
        q.solve()
        W = extract_solution(q,w)
        return (True, W)
    except MIPSolverException:
        return (False, None)


def fully_carried_solution(N):
    """
    Look for a face vector w with N*w = 0 and with all entries
    positive, minimizing the sum of its entries.  If one exists,
    returns (True, w).
    """
    num_faces = N.dimensions()[1]
    # q = MixedIntegerLinearProgram( maximization = False, solver = 'Gurobi' ) ### Grrr.
    # q = MixedIntegerLinearProgram( maximization = False, solver = 'Coin' ) ### Why!!!
    # q = MixedIntegerLinearProgram( maximization = False, solver = 'GLPK' )
    q = MixedIntegerLinearProgram( maximization = False, solver = 'PPL' )
    w = q.new_variable(real = True, nonnegative = True)
    # w = q.new_variable(integer = True, nonnegative = True)
    for v in N.rows():
        q.add_constraint( dot_prod(v, w) == 0 )
    for i in range(num_faces):
        q.add_constraint( w[i] >= 1 )
    S = sum( w[i] for i in range(num_faces) ) # Euler char
    # S = q.sum( w[i] for i in range(num_faces) ) # Euler char
    q.set_objective( S )
    try:
        q.solve()
        W = extract_solution(q,w)
        out = (True, W)
    except MIPSolverException:
        out = (False, None)
    return out


@liberal
def is_layered(tri, angle):
    """
    Given a tri, angle, decide if it is layered.
    """
    N = edge_equation_matrix_taut(tri, angle)
    N = Matrix(N)
    layered, _ = fully_carried_solution(N)
    return layered


# TODO - check this against the homological dim of the taut cone... and then delete.  :)
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
        out = 'L'
    elif non_triv:
        out = 'M'
    else:
        assert farkas
        out = 'N'
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
        assert non_triv or not full ### full => non_triv
        assert alex or not full ### full => fibered => alex is monic
        assert hyper or not full ### full => fibered => hyper is monic
        assert alex or not hyper ### hyper is monic => alex is monic
    except AssertError:
        print 'contradiction in maths'
        raise
    if full:
        pass
    elif non_triv:
        print ' non-triv sol (but not full)',
        print alex, hyper
    elif not alex or not hyper:
        print ' no sol',
        print alex, hyper
    return None

# Code to compute the dimension of the taut cone as projected into
# homology.

# (co)homology


def matrix_transpose(M):
    return map(lambda *row: list(row), *M)


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
            perm = tet.faceMapping(2,i)
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
            perm = tri.edgeMapping(i) # == tri.faceMapping(1,i)
            row[edge_index] += perm.sign()
        matrix.append(row)
    return matrix_transpose(matrix)


# extracting the taut cone


def elem_vector(i, dim):
    vec = [0]*dim
    vec[i] = 1
    return vec


def taut_rays(N):
    # get the extreme rays of the taut cone
    dim = N.dimensions()[1]
    elem_ieqs = [[0] + elem_vector(i, dim) for i in range(dim)]
    N_rows = [v.list() for v in N.rows()]
    N_eqns = [[0] + v for v in N_rows]
    P = Polyhedron(ieqs = elem_ieqs, eqns = N_eqns)
    return [ray.vector() for ray in P.rays()]


# the function


@liberal
def taut_cone_homological_dim(tri, angle):
    # find the dimension of the projection of the taut cone into
    # homology

    # boundaries of tets
    bdys = zeroth_coboundary(tri)
    bdys = matrix_transpose(bdys)

    N = edge_equation_matrix_taut(tri, angle)
    N = Matrix(N)
    rays = taut_rays(N)
    # but these are all "upwards", so we need to fix the
    # co-orientations
    coorient = is_transverse_taut(tri, angle, return_type = "face_coorientations")
    rays = [[int(a*b) for a, b in zip(coorient, ray)] for ray in rays]

    # now work in the space of two-chains
    Rays = IntegerLattice(rays + bdys)
    Bdys = Matrix(bdys)
    Cobs = Bdys.transpose()
    Anns = Cobs.kernel()
    return Rays.intersection(Anns).dimension()


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
