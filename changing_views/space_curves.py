"""Visual events for a space curve C = {F1 = F2 = 0} in P^3 (complete intersection).

F1, F2 are homogeneous polynomials in the same ring with four variables over a
field K. Points of C are handled in the chart x0 = 1.
"""
from sage.all import PolynomialRing, matrix, vector

from .binary_forms import common_roots_conditions
from .grassmannian import (chart_closure, lines_through_point, plucker_ring,
                           restrict_to_lines)


def _gradients(F1, F2, x):
    """2x4 matrix of the gradients of F1, F2 evaluated at the point x."""
    R = F1.parent()
    return matrix([[F.derivative(xi)(*x) for xi in R.gens()] for F in (F1, F2)])


def tangential_surface(F1, F2):
    """The union of the tangent lines of C, as a surface in P^3 (coordinates y)."""
    K = F1.base_ring()
    S = PolynomialRing(K, ['x1', 'x2', 'x3', 'y0', 'y1', 'y2', 'y3'])
    x = vector(S, [1] + list(S.gens()[:3]))
    y = vector(S, S.gens()[3:])
    I = S.ideal([F1(*x), F2(*x)] + list(_gradients(F1, F2, x) * y))
    I = I.elimination_ideal(S.gens()[:3])
    return I.change_ring(PolynomialRing(K, ['y0', 'y1', 'y2', 'y3']))


def trisecant_curve(F1, F2):
    """Lines meeting C in three points, as a curve in Gr(1, P^3)."""
    a1, a2 = restrict_to_lines(F1), restrict_to_lines(F2)
    A = a1[0].parent()
    return chart_closure(A.ideal(common_roots_conditions(a1, a2, 3)))


def edge_curve(F1, F2):
    """Lines through two points of C whose tangent lines meet (equivalently,
    lines lying in a plane tangent to C at two points), as a curve in Gr(1, P^3)."""
    K = F1.base_ring()
    G = plucker_ring(K)
    S = PolynomialRing(K, ['x1', 'x2', 'x3', 'y0', 'y1', 'y2', 'y3'] + list(G.variable_names()))
    x = vector(S, [1] + list(S.gens()[:3]))
    y = vector(S, S.gens()[3:7])
    p = S.gens()[7:]
    tangents = _gradients(F1, F2, x).stack(_gradients(F1, F2, y))
    minors = matrix(S, [x, y]).minors(2)
    I = S.ideal([F1(*x), F2(*x), F1(*y), F2(*y), tangents.det()]
                + [p[k] - minors[k] for k in range(6)])
    I = I.elimination_ideal(S.gens()[:7])
    return I.change_ring(G)


def chow_hypersurface(F1, F2):
    """Lines meeting C, as a hypersurface section of Gr(1, P^3). Its degree is
    twice the degree of the Chow form (which equals deg C)."""
    K = F1.base_ring()
    G = plucker_ring(K)
    S = PolynomialRing(K, ['x1', 'x2', 'x3'] + list(G.variable_names()))
    x = vector(S, [1] + list(S.gens()[:3]))
    I = S.ideal([F1(*x), F2(*x)] + lines_through_point(S.gens()[3:], x))
    I = I.elimination_ideal(S.gens()[:3])
    return I.change_ring(G)


def secant_congruence(chow):
    """Lines meeting C twice: the singular locus of the Chow hypersurface."""
    G = chow.ring()
    J = matrix([[g.derivative(v) for v in G.gens()] for g in chow.gens()])
    return chow + G.ideal(J.minors(2))
