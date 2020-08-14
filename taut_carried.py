#
# taut_carried.py
#

# functions for working with surfaces carried by the two-skeleton of a 
# transverse taut triangulation.

import regina # needed inside of imported files
import taut
import transverse_taut

def boundary_cycles_from_surface(tri, angle, face_coorientations, surface):
    """ Takes a carried surface. For each cusp of tri, look at the boundary curve
    of the surface on the boundary torus for that cusp. Push it up slightly, record 
    which faces of tri it goes through."""

    ### set up output vectors
    out = []
    for vertex in tri.faces(0):   ## 0 is the dimension of the face, so this is cusps
        out.append([0] * tri.countFaces(2))

    for f in tri.faces(2):
        for vert_num in range(3):
            cusp_index = f.face(0, vert_num).index()

            sgn = face_coorientations[f.index()]
            trailing_edge = f.face(1,(vert_num + sgn)%3) 
            leading_edge = f.face(1,(vert_num - sgn)%3) 
            ## these are with respect to right hand rule on surface as viewed from above
            ## the faces above us on the trailing edge get -1s, 
            ## the faces above us on the leading edge get +1s.

            

