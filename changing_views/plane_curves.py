"""Flexes, bitangents and the dual of a plane curve C = {f = 0} in P^2.

f is a form in three variables x0, x1, x2 over QQ. Lines u0*x0 + u1*x1 + u2*x2 = 0
are points of the dual plane with coordinates u0, u1, u2. Computations use a
random change of coordinates, so that nothing special happens at infinity.
"""
from sage.all import (QQ, ComplexField, PolynomialRing, RealField, ZZ, matrix,
                      random_matrix, vector)

from .binary_forms import coincident_root_locus, pullback


def random_plane_curve(R, d):
    """A random form of degree d in the polynomial ring R."""
    K = R.base_ring()
    return sum(K.random_element() * m for m in (sum(R.gens())**d).monomials())


def trott_quartic(R):
    """Trott's quartic: four ovals and all 28 bitangents real."""
    x0, x1, x2 = R.gens()
    return 144 * (x1**4 + x2**4) - 225 * (x1**2 + x2**2) * x0**2 + 350 * x1**2 * x2**2 + 81 * x0**4


def _peanut(R, c):
    """(x^2 - 1)^2 + y^4 + 4y^2 = c in the chart x0 = 1."""
    x0, x1, x2 = R.gens()
    return (x1**2 - x0**2)**2 + x2**4 + 4 * x2**2 * x0**2 - c * x0**4


def _nested(R):
    """Two nested ellipses, perturbed into a smooth curve with two nested ovals."""
    x0, x1, x2 = R.gens()
    return (x1**2 + 2 * x2**2 - x0**2) * (4 * x1**2 + 9 * x2**2 - 36 * x0**2) + x0**4


# smooth real quartics with different real configurations; by Klein's formula
# 2 * (isolated real bitangents) + (real flexes) = 8 for each of them
EXAMPLES = {
    'trott': trott_quartic,                                   # 4 ovals, 28 real bitangents
    'two ovals': lambda R: _peanut(R, QQ(9) / 10),             # 4 + 4 isolated bitangents
    'peanut': lambda R: _peanut(R, QQ(3) / 2),                 # 4 real flexes
    'nested ovals': _nested,                                   # 4 isolated bitangents
}


def example_curve(R, name):
    return EXAMPLES[name](R)


def _generic_matrix(n=3):
    while True:
        M = random_matrix(ZZ, n, n, x=-5, y=6)
        if M.det() != 0:
            return M


def dual_ring(K):
    return PolynomialRing(K, 'u0,u1,u2')


def hessian_curve(f):
    R = f.parent()
    return matrix(R, 3, 3, lambda i, j: f.derivative(R.gen(i), R.gen(j))).det()


def flex_points(f):
    """Ideal of the inflection points: C intersected with its Hessian curve."""
    return f.parent().ideal([f, hessian_curve(f)])


def _lines_in_chart(f, pattern, M):
    """Lines with the given contact, computed in the chart of lines not passing
    through the point M*(1, 0, 0). Isolated lines through that point are lost."""
    R = f.parent()
    K = R.base_ring()
    U = dual_ring(K)
    # lines u' = (1, u1, u2) in the coordinates x = M x', spanned by
    # x' = (-u1 - u2*t, 1, t)
    A = PolynomialRing(K, 'u1,u2')
    u1, u2 = A.gens()
    t = PolynomialRing(A, 't').gen()
    a = f(*(M * vector([-u1 - u2 * t, 1, t]))).list()
    a += [0] * (f.degree() + 1 - len(a))
    I = A.ideal(pullback(coincident_root_locus(f.degree(), pattern), a))
    # projective closure in the dual plane, then back to the old coordinates
    affine = U.ideal([g(*U.gens()[1:]) for g in I.gens()])
    closure = [g.homogenize(U.gen(0)) for g in affine.groebner_basis()]
    back = M.transpose() * vector(U.gens())  # u' = M^T u
    return U.ideal([g(*back) for g in closure])


def lines_with_contact(f, pattern):
    """Ideal in the dual plane of the lines meeting C with the given contact
    pattern: '2' gives the dual curve, '3' the flex lines, '2+2' the bitangents.

    Computed in two random charts, so that no line is lost: each chart misses
    only isolated lines through one random point."""
    from .grassmannian import dimension_and_degree
    I1 = _lines_in_chart(f, pattern, _generic_matrix())
    I2 = _lines_in_chart(f, pattern, _generic_matrix())
    if dimension_and_degree(I1) == dimension_and_degree(I2):
        return I1
    return I1.intersection(I2)


def dual_curve(f):
    return lines_with_contact(f, '2')


def flex_lines(f):
    return lines_with_contact(f, '3')


def bitangents(f):
    return lines_with_contact(f, '2+2')


def real_points(I, prec=100):
    """Real points of a zero-dimensional radical homogeneous ideal in three
    variables, as vectors over RealField(prec) normalized to unit length.

    In random coordinates the lex Groebner basis has the shape {a - g(b), h(b)}
    with deg h = deg I. The count is exact (real root isolation of h); the
    coordinates are computed with ample working precision."""
    from .grassmannian import dimension_and_degree
    S = I.ring()
    A = PolynomialRing(S.base_ring(), 'a,b', order='lex')
    a, b = A.gens()
    n = dimension_and_degree(I)[1]
    while True:  # random coordinates until the shape lemma applies
        M = _generic_matrix()
        chart = M * vector(A, [1, a, b])
        gb = A.ideal([g(*chart) for g in I.gens()]).groebner_basis()
        if (len(gb) == 2 and gb[0].degree(a) == 1 and gb[0].lc() == 1
                and gb[0].coefficient({a: 1}) == 1 and gb[1].degree() == n):
            break
    h = gb[1].univariate_polynomial()
    g = (a - gb[0]).univariate_polynomial(h.parent())
    work = RealField(max(10 * prec, 1000))
    RF = RealField(prec)
    points = []
    for r in h.roots(work, multiplicities=False):
        p = M.change_ring(work) * vector(work, [1, g(r), r])
        points.append((p / p.norm()).change_ring(RF))
    return points


def tangent_line(f, p):
    """The tangent line (as a vector u) at a point p of C."""
    R = f.parent()
    return vector([f.derivative(x)(*p) for x in R.gens()])


def contact_points(f, u, prec=100):
    """The intersection points of C with the line u, as (point, is_real) pairs."""
    RF, CF = RealField(prec), ComplexField(prec)
    u = vector(RF, u)
    u = u / u.norm()
    # orthonormal basis P, Q of the line, well conditioned for any u
    e = vector(RF, [0, 0, 0])
    e[min(range(3), key=lambda i: abs(u[i]))] = 1
    P = u.cross_product(e)
    P = P / P.norm()
    Q = u.cross_product(P)
    t = PolynomialRing(RF, 't').gen()
    g = f.change_ring(RF)(*(P + t * Q))
    points = []
    for r in g.roots(CF, multiplicities=False):
        # a real double root splits by about sqrt(precision)
        real = abs(r.imag()) < 1e-8 * (1 + abs(r))
        points.append((P + r.real() * Q if real else P + r * Q, real))
    return points


def is_isolated_bitangent(f, u, prec=100):
    """True if the real bitangent u touches C at two non-real points."""
    return not any(real for _, real in contact_points(f, u, prec))


def summary(f):
    """Complex and real counts of flexes and bitangents, with the checks given by
    the Plücker formulas (for smooth curves) and Klein's formula
        d + 2 * (isolated real bitangents) + (real flexes) = d(d - 1)."""
    from .grassmannian import dimension_and_degree
    d = f.degree()
    flexes, bitang = flex_points(f), bitangents(f)
    real_flexes = real_points(flexes)
    real_bitangents = real_points(bitang)
    isolated = [u for u in real_bitangents if is_isolated_bitangent(f, u)]
    s = {
        'degree': d,
        'class': dimension_and_degree(dual_curve(f))[1],
        'flexes': dimension_and_degree(flexes)[1],
        'bitangents': dimension_and_degree(bitang)[1],
        'real flexes': len(real_flexes),
        'real bitangents': len(real_bitangents),
        'isolated real bitangents': len(isolated),
    }
    s['Plücker'] = (s['class'] == d * (d - 1) and s['flexes'] == 3 * d * (d - 2)
                    and s['bitangents'] == d * (d - 2) * (d - 3) * (d + 3) // 2)
    s['Klein'] = d + 2 * s['isolated real bitangents'] + s['real flexes'] == s['class']
    return s
