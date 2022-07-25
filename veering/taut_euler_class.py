#
# taut_euler_class.py
#

from sage.matrix.constructor import Matrix
from sage.modules.free_module_element import vector
from sage.arith.misc import gcd
from sage.arith.functions import lcm

from .file_io import parse_data_file, write_data_file
from .taut import liberal, isosig_to_tri_angle
from .transverse_taut import is_transverse_taut


#
# Goal - given a transverse taut triangulation, decide if the
# associated "absolute" euler class is torsion or not.  If it is
# torsion, determine its order.
#


# Contents and overview:


# 1. References.
#
# 2. Background.
#
# 3. Helper functions.
#
# 4. Truncate.  We build the correct "truncated" cell structure \calT'
# from (M, \calT) and give generators for the cochain groups
# C^k(\calT', \ZZ) (k = 1, 2).
#
# 5. Representative.  We find a two-cocycle E \in Z^2(\calT', \ZZ)
# that represents E(\calT) \in H^2(M, \ZZ).
#
# 6. Coboundary.  We find the matrix for the coboundary operator
# \delta^1.
#
# 7. Linear algebra.  We solve the linear problem to decide if E is a
# coboundary - that is, if E lies in B^2(\calT', \ZZ) - that is, if E
# is in the image of \delta^1.
#
# 8. Remarks.
#
# 9. Calling code
#


# 1. References.
#
# Culler, Dunfield - Orderability and Dehn filling
# Ghys - Groups acting on the circle
# Thurston - A norm for the homology of three-manifolds
# Candel, Conlon - Foliations, chapter four


# 2. Background:

# Suppose that (M, \calT) is a transverse taut triangulation.  Then
# \calT^{2} is the "horizontal branched surface".  This caries various
# laminations, which extend to foliations on M.  All of these have the
# same Euler class, which we will denote E(\calT) \in H^2(M, \ZZ).
# Suppose that \calF is a carried foliation and let UT\calF be the
# unit tangent bundle over \calF.  The Euler class E vanishes exactly
# when UT\calF has a section; that is, when the unit tangent bundle is
# trivialisable.


# Recall:
# Suppose that X is an F-bundle over B.  We have
#
#      i
# F -------> X <--.
#            |    |
#            |    |
#           p|    |s
#            |    |
#            v    |
#            B ---'
#
# So s \from B \to X is a \emph{section} if p \circ s = Id_B


# 3. Helper functions


def diagonal(D):
    return [D[i][i] for i in range(min(D.dimensions()))]


# 4. Truncate.

# Suppose that M is a connected, cusped, oriented three-manifold.  Let
# C = C(M) \geq 1 be the number of cusps of M.  Suppose that \calT is a
# transverse taut ideal triangulation of M.  Let T = T(\calT) \geq 1
# be the number of tetrahedra of \calT.

# We use Regina to number and orient the edges \{e_i\}_{i = 0}^{T-1},
# the faces \{f_i\}_{i = 0}^{2T-1}, and the tetrahedra \{t_i\}_{i =
# 0}^{T-1} of \calT.  We call all of these \emph{ideal} cells.  Note
# that the transverse structure also gives us co-orientations of the
# e_i and the f_i, called "upwards"

# We remove a small open neighbourbood of all ideal vertices of all
# model tetrahedra.  This gives the \emph{truncated} cell structure
# \calT'.  The remains of the ideal cells are called \emph{truncated}
# cells; we abuse and reuse the notations e_i and f_i for these.  The
# truncated cells inherit orientations and co-orientations.  The new
# cells are called \emph{peripheral} cells.  We number these as
# follows:

# e_{ij} is the peripheral edge cutting vertex v_j off of ideal face f_i
# f_{ij} is the peripheral face cutting vertex v_j off of ideal tetrahedron t_i

# Note that every truncated face is combinatorially a hexagon.  The
# boundary of this hexagon contains three truncated edges alternating
# with three peripheral edges.  We orient each peripheral edge e_{ij}
# so that the orientation of e_{ij} agrees with the orientation
# induced by \bdy f_i.  We orient each peripheral face f_{ij}
# anti-clockwise, as viewed from infinity (that is, from outside of
# M).  Also, we equip e_{ij} and f_{ij} with co-orientations pointing
# out of M, called "outward".

#             e_{i0}
#              ---
#             /   \
#       e_2  /     \  e_1
#           /       \
#          /   f_i   \
#          \         /
#   e_{i1}  ---------  e_{i2}
#              e_0

# For an edge e or a face f we use e^* and f^* to denote the dual in
# C^1(\calT', \ZZ) or C^2(\calT', \ZZ).  Thus \{e^*_i\} \cup
# \{e^*_{ij}\} generates C^1(\calT', \ZZ) while \{f^*_i\} \cup
# \{f^*_{ij}\} generates C^2(\calT', \ZZ).

# For more pictures, see
# /Veering_code/NotesPictures/euler_notes_from_nathan.jpg


# 5. Representative

# We now construct a two-cocycle E \in Z^2(\calT', \ZZ).  For every
# peripheral face f we take

# E(f) = 0.

# \begin{remark}
# To see that this is correct, let \calF be any foliation of M,
# transverse to the boundary.  Suppose that f is the given peripheral
# triangle.  We have a section of the restriction of UT\calF to \bdy
# f; namely the outward field.  This extends over f to give a section
# of UT\calF restricted to f.  So there is no obstruction to the
# extension.  See below for a more precise discussion in terms of
# "Poincar\'e-Hopf index".
# \end{remark}

# Now suppose that f is a truncated face.  Suppose that e_0, e_1, e_2
# are its three truncated edges.  Recall that these are all oriented.
# Let AC(f) be the number of the edges e_0, e_1, e_2 that are
# oriented anti-clockwise (that is, agree with their induced
# orientation coming from f).  We take

# E(f) = AC(f) - 2

# If we flip the transverse direction: AC(f') = 3 - AC(f),
# so E(f') = AC(f') - 2 = 1 - AC(f) = 2 - AC(f) - 1 = -E(f) - 1

# \begin{remark}
# Here is one way to remember (and explain!) this rule.  Suppose that
# f is the given truncated face.  Suppose that s is a section of UTf |
# \bdy f.  Then index(s) is the total rotation of s with respect to
# the tangent field, _plus_ one.  This can be rephrased in terms of
# the index of tangent vector fields extending s over all of f.

# Our choices of orientations of edges determine a section of UTf |
# \bdy f.  Since all of the boundary edges e_{ij} of f are oriented
# the same way, we choose a standard framing there; Nathan tells us to
# just use the outward pointing section on all of the e_{ij}.  Our
# choice of section on e_0 (say) has to (a) depend only on the
# orientation of e_0 and (b) has to be outward at the endpoints of
# e_0.  The simplest choice is the section that rotates by +\pi with
# respect to the tangent along \bdy f_i, as we move forward along e_0.
# So s points _back_ at the beginning of e_0, points _right_ in the
# middle of e_0, and points _forwards_ at the end of e_0.  The total
# rotation of the resulting field (with respect to the tangent field)
# is AC(f) - 3.  Thus E(f) = AC(f) - 2 is the index.  You can check
# this works by drawing the four possible pictures and computing the index
# of any extension of s over f.
# \end{remark}

# Claim: \delta^2 E = 0.

# That is, E is a cocycle.

# Proof of claim: Fix a truncated tetrahedron t and fix some oriention
# of its truncated edges.  A direct calculation shows that

# \delta E (t) = E \bdy t = 0.

# Likewise, a direct computation shows that switching the orientation
# of a single edge leaves E \bdy t unchanged. QED.

### It would be nice to have a less computational proof!


def euler_cocycle(tri, angle):
    """
    Given a regina triangulation "tri", with oriented edges, and a
    transverse taut angle structure "angle", returns the associated
    two-cocycle E representing the Euler class E(tri).
    """
    assert is_transverse_taut(tri, angle)
    face_coorientations = is_transverse_taut(tri, angle, return_type = "face_coorientations")
    # E will be a _row_ vector, because it eats column vectors.
    E = []
    # First deal with the truncated faces
    for face in tri.faces(2):  # 2 = dimension
        # First we compute the number of Regina oriented edges that agree with the Regina orientation on face
        AC = 0
        for i in range(3):
            perm = face.faceMapping(1, i)
            # print perm[0], perm[1]
            if perm[1] == ((perm[0] + 1) % 3):  # the edge and face orientations agree so,
                AC = AC + 1
        # print "AC", AC
        # Now we condition on whether or not Regina and angle agree on the (co-)orientation of the face.
        if face_coorientations[face.index()] == 1:
            E.append(AC - 2)
        else:
            E.append(1 - AC)
    # Now deal with the peripheral faces
    for tet in tri.tetrahedra():
        for j in range(4):
            E.append(0)
    return E


# 6. Coboundary

# Suppose that e is a truncated edge.  Let LF be the set of truncated
# faces to the left of e and let RF be the set of faces to the right.  Then

# \delta e^* = \sum_{f \in LF} f^* - \sum_{f \in RF} f^*.

# Suppose that e is a peripheral edge.  So there is a unique truncated
# face f meeting e.  Note that f is to the left of e.  There are
# also a pair of boundary faces meeting e: say f' _above_ e and f''
# _below_ e.  Then

# \delta e^* = f^* + (f')^* - (f'')^*.


def coboundary(tri, angle):
    """
    Given a triangulation "tri" (T), with oriented edges, and a
    transverse taut angle structure "angle", returns the co-boundary
    operator delta^1 \from C^1(T', ZZ) \to C^2(T', ZZ), as a matrix,
    for the truncated triangulation T'.  Note that, strictly speaking,
    we don't need to use "angle" for this, but we use it to determine
    orientation on faces for the Euler class, so we might as well use
    it again here.
    """
    # \delta^1 takes row vectors (functions on edges) and spits out
    # row vectors (functions on faces).  So, if c is a one-cochain
    # then c \cdot \delta is a two-cochain.
    delta = []
    assert is_transverse_taut(tri, angle)
    tet_vert_coorientations = is_transverse_taut(tri, angle, return_type = "tet_vert_coorientations")
    face_coorientations = is_transverse_taut(tri, angle, return_type = "face_coorientations")

    for edge in tri.edges():
        # A row for every truncated edge
        row = []
        for face in tri.triangles():
            # A row entry for every truncated face
            count = 0
            for i in range(3):
                if face.edge(i) == edge:
                    perm = face.faceMapping(1, i)
                    if perm[1] == ((perm[0] + 1) % 3):
                        # the edge and face orientations agree so,
                        count += 1
                    else:
                        count -= 1
            row.append(count * face_coorientations[face.index()])
            # +1 if face is to the left of the edge, -1 if face is to
            # the right of the edge, using Regina's edge orientation
            # when viewed from above (using the transverse taut notion
            # of up)

            #        ,'|
            #      ,'  |
            #    ,'    |
            #  ,'  CCW |  gets a +1
            #  `.      ^
            #    `.    |
            #      `.  |
            #        `.|

        for tet in tri.simplices():
            for i in range(4):
                row.append(0)
        delta.append(row)

    for face in tri.triangles():
        face_embeddings = []
        for j in range(2):
            face_embeddings.append( face.embedding(j) )

        for i in range(3):  # vertices of the face
            # A row for every peripheral edge
            row = []

            for face2 in tri.triangles():
                # A row entry for every truncated face
                if face2 == face:
                    row.append(1)
                else:
                    row.append(0)

            for tet in tri.simplices():
                for k in range(4):
                    # A row entry for every peripheral face
                    count = 0
                    for j in range(2):
                        if (tet == face_embeddings[j].simplex()) and (face_embeddings[j].vertices()[i] == k):
                            # the tetrahedron is on the jth side of the
                            # face and the ith vertex of face is the kth
                            # vertex of tet
                            face_num_in_tet = face_embeddings[j].vertices()[3]
                            count -= tet_vert_coorientations[tet.index()][face_num_in_tet]
                            # tet_vert_coorientations is +1 if
                            # coorientation on face points out of the
                            # tetrahedron, and we want count += 1 if
                            # the peripheral face is above the
                            # peripheral edge
                    row.append(count)
            delta.append(row)
    return delta


# 7. Linear algebra

# We ask: is there a one-cocycle C \in C^1(\calT', \ZZ) so that
# \delta C = E?  If so, then [E] = E(\calT) is zero in H^2, as
# desired.

# This is a linear algebra problem, so can be solved by, say, sage.


def order_of_euler_class(delta, E):
    """
    Given the coboundary operator delta and an Euler two-cocycle E,
    returns k if [E] is k--torsion.  By convention, returns zero if
    [E] is non-torsion.  Note that the trivial element is 1--torsion.
    """

    delta = Matrix(delta)
    E = vector(E)

    # Note that E is a coboundary if there is a one-cocycle C solving
    #
    # E = C*delta
    #
    # We can find C (if it exists at all) using Smith normal form.

    D, U, V = delta.smith_form()
    assert D == U*delta*V

    # So we are trying to solve
    #
    # C*delta = C*U.inverse()*D*V.inverse() = E
    #
    # for a one-cochain C.  Multiply by V to get
    #
    # C*delta*V = C*U.inverse()*D = E*V
    #
    # Now set
    #
    # B = C*U.inverse(), and so B*U = C
    #
    # and rewrite to get
    #
    # B*U*delta*V = B*D = E*V
    #
    # So define E' by:

    Ep = E*V

    # Finally we attempt to solve B * D = Ep.  Note that D is
    # diagonal: so if we can solve all of the equations

    # B[i] * D[i][i] == Ep[i]

    # with B[i] integers, then [E] = 0 in cohomology.

    diag = diagonal(D)

    if any( (diag[i] == 0 and Ep[i] != 0) for i in range(len(Ep)) ):
        return 0

    # All zeros are at the end in Smith normal form.  Since we've
    # passed the above we can now remove them.

    first_zero = diag.index(0)
    diag = diag[:first_zero]
    Ep = Ep[:first_zero]

    # Since diag[i] is (now) never zero we can divide to get the
    # fractions Ep[i]/diag[i] and then find the scaling that makes
    # them simultaneously integral.

    denoms = [ diag[i] / gcd(Ep[i], diag[i]) for i in range(len(Ep)) ]
    return lcm(denoms)


# 8. Remarks

# a) Here is a nice trick that proves [E] = 0 in some cases.  Suppose
# that \gamma is an oriented path in \bdy M.  Suppose that \gamma is
# transverse to the one-skeleton of \calT'.  We form a one-cocycle
# D_\gamma by adding up the boundary edges that \gamma crosses, with
# sign.  The sign is positive if \gamma crosses from below to above,
# and negative otherwise.  Note that \delta D_\gamma vanishes on all
# boundary faces.

# b) Marc Lackenby says that we should take the paths that go up
# through the centres of tetrahedra and take the Poincare dual.  BUT I
# think this is not what we want... Marc is thinking of the relative
# Euler class as discussed on page 390 of his paper "Taut ideal
# triangulations of three-manifolds".  The relative Euler class lives
# in H^2(M, \bdy M), so is Poincare dual to an element of H_1(M),
# represented by a collection of loops.

# c) [2019-03-31] It seems that, for transverse veering triangulations
# in the 16 census, the Euler class is always zero or two-torsion.
# Note that there are manifolds M in the census where H^2(M, \ZZ) has
# positive rank...  What about odd torsion?

# Question: If the veering triangulation is edge-orientable, does the
# Euler class vanish?

# Answer: Yes.  Here is a version of a discussion with Nathan
# [2020-04-03] - he says the following:

# Suppose that F is a foliation carried by the horizontal branched
# surface.  Let UTF be the unit tangent bundle to F.  We think of
# e(UTF) as being the obstruction to UTF having a section.  Let G be
# the foliation carried by the upper (aka green) branched surface.  If
# G is transversely orientable (aka edge-orientability of the veering
# triangulation) then G \cap F gives the desired section, and e(UTF) =
# 0.  Note that G \cap F gives, for every point, a pair of points in
# the unit tangent circle.  So let PUTF be the projective unit tangent
# bundle to F.  This definitely has a section, so e(PUTF) = 0.  Now,
# the bundle UTF is a double cover of the bundle PUTF.

# Claim: The euler class is multiplicative with respect to covers (in
# both senses).

# With the claim in hand, we have

# 2 * e(UTF) = e(PUTF) = 0

# We deduce that e(UTF) is either zero or two-torsion.


# 9. Calling code


@liberal
def order_of_euler_class_wrapper(tri, angle):
    """
    Returns the order of the euler class.
    """
    return order_of_euler_class(coboundary(tri, angle), euler_cocycle(tri, angle))


def compute_order_of_euler_classes(file_in, number=None, file_out=None):
    data_in = parse_data_file(file_in)
    data_in = [line.split(" ") for line in data_in]
    if number != None:
        data_in = data_in[:number]
    data_out = []
    evil = []
    for i, line in enumerate(data_in):
        if i % 50 == 0:
            print( ((1.0*i)/(1.0*len(data_in)), len(data_out)) )
        sig = line[0]
        tri, angle = isosig_to_tri_angle(sig)
        # angle = [int(letter) for letter in angle_s]
        curr_euler = order_of_euler_class(coboundary(tri, angle), euler_cocycle(tri, angle))
        if curr_euler == "non-torsion":
            evil.append(sig)
            print(sig + " has non-torsion Euler class!!!!")
        elif curr_euler == 1:  # order is one so [E] = 0.  Boring.
            pass
        else:
            line_out = [sig, str(curr_euler)]
            line_out.extend(line[1:])
            data_out.append(line_out)
    if file_out != None:
        write_data_file(data_out, file_out)
    print( ("list of evil:", evil) )
    return data_out
