"""Checks against classical formulas. Run with

    sage -python tests/test_events.py

(or with pytest, if installed in the Sage environment).
"""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from sage.all import GF, PolynomialRing, QQ, set_random_seed

from changing_views.binary_forms import coincident_root_locus
from changing_views.grassmannian import dimension_and_degree, order_and_class
from changing_views.space_curves import (chow_hypersurface, secant_congruence,
                                        trisecant_curve)
from changing_views.surfaces import (congruence, contour, event_curve,
                                    hurwitz_threefold, parabolic_curve,
                                    random_surface)

K = GF(32003)
R = PolynomialRing(K, 'x0,x1,x2,x3')


def test_discriminant():
    # classical discriminant of a quadratic, and the twisted cubic of cubes
    a0, a1, a2 = coincident_root_locus(2, '2').ring().gens()
    assert coincident_root_locus(2, '2').gens() == [a1**2 - 4 * a0 * a2]
    I = coincident_root_locus(3, '3')
    assert I.is_prime() and dimension_and_degree(I) == (1, 3)


def test_quartic_surface():
    set_random_seed(0)
    d = 4
    F = random_surface(R, d)
    assert dimension_and_degree(contour(F, [1, 2, 3, 5])) == (1, d * (d - 1))
    assert dimension_and_degree(hurwitz_threefold(F)) == (3, 2 * d * (d - 1))
    assert order_and_class(congruence(F, 'bitangent')) == (
        d * (d - 1) * (d - 2) * (d - 3) // 2, d * (d - 2) * (d - 3) * (d + 3) // 2)
    assert order_and_class(congruence(F, 'principal')) == (
        d * (d - 1) * (d - 2), 3 * d * (d - 2))
    # Table 2 of Kohn-Sturmfels-Trager, arXiv:1707.01877
    assert dimension_and_degree(event_curve(F, 'flecnodal')) == (1, 80)


def test_cubic_surface():
    set_random_seed(0)
    # Table 2 of Kohn-Sturmfels-Trager
    assert dimension_and_degree(parabolic_curve(random_surface(R, 3))) == (1, 30)


def test_space_sextic():
    set_random_seed(0)
    F1, F2 = random_surface(R, 2), random_surface(R, 3)
    n, g = 6, 4
    assert dimension_and_degree(trisecant_curve(F1, F2)) == (
        1, (n - 2) * ((n - 1) * (n - 3) - 3 * g) // 3)
    chow = chow_hypersurface(F1, F2)
    assert dimension_and_degree(chow) == (3, 2 * n)
    assert order_and_class(secant_congruence(chow)) == (
        (n - 1) * (n - 2) // 2 - g, n * (n - 1) // 2)


if __name__ == '__main__':
    for name, f in list(globals().items()):
        if name.startswith('test_'):
            f()
            print(name, 'ok')
