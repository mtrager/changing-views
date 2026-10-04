"""Lines in P^3: Plücker coordinates, restriction to lines, and invariants.

Lines are represented in the Plücker embedding Gr(1, P^3) in P^5, with
coordinates p01, p02, p03, p12, p13, p23 where p_ij = x_i y_j - x_j y_i for a
line through points x and y. Computations happen in the affine chart p01 = 1,
where the line is spanned by P = (1, 0, -p12, -p13) and Q = (0, 1, p02, p03).
"""
from sage.all import PolynomialRing, factorial, matrix, vector

from .binary_forms import coincident_root_locus, pullback

PLUCKER = ('p01', 'p02', 'p03', 'p12', 'p13', 'p23')
CHART = ('p02', 'p03', 'p12', 'p13')


def plucker_ring(K):
    return PolynomialRing(K, PLUCKER)


def chart_ring(K):
    return PolynomialRing(K, CHART)


def plucker_relation(G):
    p01, p02, p03, p12, p13, p23 = G.gens()
    return p01 * p23 - p02 * p13 + p03 * p12


def chart_line(K):
    """The points P, Q spanning the line with chart coordinates (p02, p03, p12, p13)."""
    p02, p03, p12, p13 = chart_ring(K).gens()
    return vector([1, 0, -p12, -p13]), vector([0, 1, p02, p03])


def restrict_to_lines(F):
    """Coefficients [a0, ..., ad] of the binary form F(P + t*Q), as polynomials
    in the chart coordinates. F is a form of degree d in four variables."""
    P, Q = chart_line(F.base_ring())
    t = PolynomialRing(P.base_ring(), 't').gen()
    f = F(*(P + t * Q)).list()
    return f + [0] * (F.degree() + 1 - len(f))


def chart_closure(I):
    """Closure in Gr(1, P^3) of an ideal I in the chart coordinates.

    Returns a homogeneous ideal in the Plücker ring (containing the Plücker
    relation)."""
    G = plucker_ring(I.base_ring())
    p01, p02, p03, p12, p13, p23 = G.gens()
    affine = G.ideal([g(p02, p03, p12, p13) for g in I.gens()]
                     + [p23 - p02 * p13 + p03 * p12])
    # homogenizing a degree-compatible Groebner basis gives the projective closure
    return G.ideal([g.homogenize(p01) for g in affine.groebner_basis()])


def lines_with_contact(F, pattern):
    """Lines meeting the surface F = 0 with the given contact pattern, e.g.
    '2' (tangent), '2+2' (bitangent), '4' (four-fold contact).

    For d <= 3 the lines contained in the surface are removed (a general surface
    of degree >= 4 contains none, and the saturation is expensive)."""
    a = restrict_to_lines(F)
    A = a[0].parent()
    I = A.ideal(pullback(coincident_root_locus(len(a) - 1, pattern), a))
    if len(a) - 1 <= 3:
        I = I.saturation(A.ideal(a))[0]
    return chart_closure(I)


def lines_through_point(p, z):
    """Equations in the Plücker coordinates p = (p01, ..., p23) for the line
    to pass through the point z (an alpha-plane)."""
    p01, p02, p03, p12, p13, p23 = p
    M = matrix([[0, p23, -p13, p12], [-p23, 0, p03, -p02],
                [p13, -p03, 0, p01], [-p12, p02, -p01, 0]])
    return list(M * vector(z))


def lines_in_plane(p, w):
    """Equations in the Plücker coordinates p = (p01, ..., p23) for the line
    to lie in the plane w (a beta-plane)."""
    p01, p02, p03, p12, p13, p23 = p
    M = matrix([[0, p01, p02, p03], [-p01, 0, p12, p13],
                [-p02, -p12, 0, p23], [-p03, -p13, -p23, 0]])
    return list(M * vector(w))


def dimension_and_degree(I):
    """Projective dimension and degree of a homogeneous ideal (degree counts
    multiplicities, as given by the Hilbert polynomial)."""
    hp = I.hilbert_polynomial()
    if hp == 0:
        return -1, 0
    dim = hp.degree()
    return dim, hp.leading_coefficient() * factorial(dim)


def order_and_class(congruence):
    """Order (number of lines through a general point) and class (number of
    lines in a general plane) of a 2-dimensional family of lines."""
    G = congruence.ring()
    K = G.base_ring()
    z = [K.random_element() for _ in range(4)]
    w = [K.random_element() for _ in range(4)]
    p = G.gens()
    return (dimension_and_degree(congruence + G.ideal(lines_through_point(p, z)))[1],
            dimension_and_degree(congruence + G.ideal(lines_in_plane(p, w)))[1])
