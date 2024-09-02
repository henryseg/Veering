### Contributed by Matthias Goerner

from snappy import Manifold
try:
    from snappy.SnapPy import list_as_word
except ImportError:
    def _to_letter(g):
        if g > 0:
            return chr(64 + g)
        else:
            return chr(96 - g)
    
    def list_as_word(l, num_gens, verbose_form):
        if verbose_form:
            raise ValueError("verbose_form not supported")
        if num_gens > 26:
            raise ValueError("Only supports up to 26 generators")
        return ''.join(_to_letter(g) for g in l)

from typing import Sequence, Tuple, Optional
                     
def tet_and_face_indices_to_word(
        mfd : Manifold,
        tet_and_face_indices : Sequence[Tuple[int, int]]) -> str:

    mfd._choose_generators(False, False)
    gen_infos = [info['generators'] for info in mfd._choose_generators_info()]

    word_list = [
        g  ### replace by -g if everything explodes
        for tet_index, face_index in tet_and_face_indices
        if (g := gen_infos[tet_index][face_index]) != 0 ]

    num_gens = mfd.fundamental_group(False).num_generators()

    return list_as_word(word_list, num_gens, verbose_form = False)

def drill_tet_and_face_indices(
        mfd : Manifold,
        tet_and_face_indices : Sequence[Tuple[int, int]],
        verified : bool = False,
        bits_prec : Optional[int] = None) -> Manifold:
    """
    Drills geodesic homotopic to given curve.
    """

    return mfd.drill_word(
        tet_and_face_indices_to_word(mfd, tet_and_face_indices))

def main():
    M = Manifold("m004")
    # Drilling aBcB
    print(
        drill_tet_and_face_indices(
            M, [(1,1), (0,1), (1,3), (0,1)]).identify())

if __name__ == '__main__':
    main()
