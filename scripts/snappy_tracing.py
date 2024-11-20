
from snappy.snap.t3mlite import Mcomplex
from snappy.geometric_structure import (
    add_r13_geometry,
    add_filling_information)
from snappy.geometric_structure.geodesic.add_core_curves import add_r13_core_curves
from snappy.geometric_structure.geodesic.geodesic_start_point_info import compute_geodesic_start_point_info
from snappy.drilling.tracing import trace_geodesic
from snappy.drilling.cusps import index_geodesics_and_add_post_drill_infos
from snappy.drilling.perturb import perturb_geodesic
from snappy.len_spec.geometric_structure import _add_and_scale_cusp_cross_section
from snappy.geometric_structure.cusp_neighborhood.vertices import scale_vertices_from_horotriangles
from snappy.snap.t3mlite import simplex

from sage.all import matrix

def hyperboloid_to_barycentric(tet, r13_point):
    hyperboloid_embedding = matrix(
        [tet.R13_vertices[v]
         for v in simplex.ZeroSubsimplices ]).transpose()
    inverse = hyperboloid_embedding.inverse()
    unnormalised_bary = inverse * r13_point
    norm = sum(unnormalised_bary)
    return unnormalised_bary / norm

def my_trace_geodesic(manifold, word, perturb_amt = 0, bits_prec = 53, verified = False):
    mcomplex = Mcomplex(manifold)
    add_r13_geometry(mcomplex, manifold,
                     verified=verified, bits_prec=bits_prec)
    add_filling_information(mcomplex,manifold)
    add_r13_core_curves(mcomplex,manifold)

    _add_and_scale_cusp_cross_section(mcomplex)
    scale_vertices_from_horotriangles(mcomplex)

    g = compute_geodesic_start_point_info(mcomplex, word)

    if perturb_amt > 0:
        RF = mcomplex.Tetrahedra[0].ShapeParameters[simplex.E01].real().parent()
        rf_perturb_amt = RF(perturb_amt)
        perturb_geodesic(g, rf_perturb_amt, verified=verified)
    
    index_geodesics_and_add_post_drill_infos([g], mcomplex)

    for piece in trace_geodesic(g, verified=verified):
        yield {
            'tet' : piece.tet.Index,
            'startBary' : hyperboloid_to_barycentric(
                piece.tet, piece.endpoints[0].r13_point),
            'startCell' : piece.endpoints[0].subsimplex,
            'endBary' : hyperboloid_to_barycentric(
                piece.tet, piece.endpoints[1].r13_point),
            'endCell' : piece.endpoints[1].subsimplex }

if __name__ == '__main__':
    from snappy import Manifold
    M = Manifold("m004")
    
    print(list(my_trace_geodesic(M, 'a', perturb_amt = 0.001)))
