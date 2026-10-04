"""Visual events for a surface S = {F = 0} in P^3.

All functions take a homogeneous polynomial F in a ring with four variables
(x0, x1, x2, x3) over a field K; working over a finite field such as GF(32003)
is much faster than QQ and gives the same degrees for generic input.
"""
from sage.all import PolynomialRing, matrix, vector

from .grassmannian import plucker_ring, lines_with_contact

# contact patterns of lines with the surface
EVENT_CURVES = {'flecnodal': '4', 'cusp crossing': '3+2', 'tritangent': '2+2+2'}
CONGRUENCES = {'bitangent': '2+2', 'principal': '3'}
HURWITZ = '2'


def random_surface(R, d):
    """A random form of degree d in the polynomial ring R."""
    K = R.base_ring()
    return sum(K.random_element() * m for m in (sum(R.gens())**d).monomials())


def event_curve(F, name):
    """The flecnodal, cusp crossing or tritangent lines, as a curve in Gr(1, P^3)."""
    return lines_with_contact(F, EVENT_CURVES[name])


def congruence(F, name):
    """The bitangent or principal (inflectional) lines, as a surface in Gr(1, P^3)."""
    return lines_with_contact(F, CONGRUENCES[name])


def hurwitz_threefold(F):
    """The tangent lines, as a hypersurface section of Gr(1, P^3). Its degree is
    twice the degree of the Hurwitz form."""
    return lines_with_contact(F, HURWITZ)


def parabolic_curve(F):
    """Asymptotic tangent lines at the parabolic points, as a curve in Gr(1, P^3).

    At a parabolic point x the Hessian H(x) is singular, and its kernel y spans
    together with x the unique asymptotic line. We work in the chart x0 = 1."""
    R = F.parent()
    K = R.base_ring()
    S = PolynomialRing(K, ['x1', 'x2', 'x3', 'y0', 'y1', 'y2', 'y3'] + list(plucker_ring(K).variable_names()))
    x = vector(S, [1] + list(S.gens()[:3]))
    y = vector(S, S.gens()[3:7])
    p = S.gens()[7:]
    Fx = F(*x)
    H = matrix(S, 4, 4, lambda i, j: F.derivative(R.gen(i), R.gen(j))(*x))
    minors = matrix(S, [x, y]).minors(2)  # in the order 01, 02, 03, 12, 13, 23
    I = S.ideal([Fx, H.det()] + list(H * y) + [p[k] - minors[k] for k in range(6)])
    I = I.elimination_ideal(S.gens()[:7])
    return I.change_ring(plucker_ring(K))


def contour(F, center):
    """The image contour (apparent contour) of S seen from the point `center`,
    as a plane curve in coordinates u0, u1, u2 of the image plane."""
    R = F.parent()
    K = R.base_ring()
    proj = matrix(K, [center]).right_kernel_matrix()  # 3x4, kills the center
    S = PolynomialRing(K, list(R.variable_names()) + ['u0', 'u1', 'u2'])
    x = vector(S, S.gens()[:4])
    u = vector(S, S.gens()[4:])
    Fx = F(*x)
    polar = sum(c * Fx.derivative(xi) for c, xi in zip(center, x))
    I = S.ideal([Fx, polar] + list(proj * x - u))
    I = I.elimination_ideal(S.gens()[:4])
    return I.change_ring(PolynomialRing(K, ['u0', 'u1', 'u2']))
