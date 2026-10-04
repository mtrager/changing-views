"""Binary forms: coincident root loci and Sylvester-type rank conditions.

A binary form of degree d is written as f(t) = a0 + a1*t + ... + ad*t^d
(the dehomogenization s = 1 of sum a_i s^(d-i) t^i).
"""
from collections import Counter
from functools import lru_cache

from sage.all import QQ, PolynomialRing, matrix


def parse_pattern(pattern):
    """'3+2' -> [3, 2]. Each part is the multiplicity of one root."""
    return sorted((int(k) for k in str(pattern).split('+')), reverse=True)


@lru_cache(maxsize=None)
def coincident_root_locus(d, pattern):
    """Ideal of the binary forms of degree d having roots with the given
    multiplicities, e.g. pattern '2+2' (two double roots) or '3' (a triple root).

    This is the projective closure in P^d with coordinates a0, ..., ad, computed
    by implicitizing f = g * prod_m h_m^m, where h_m is a monic polynomial whose
    roots are the roots of multiplicity m and g is arbitrary.
    """
    groups = Counter(parse_pattern(pattern))  # multiplicity -> number of roots
    free = d - sum(m * k for m, k in groups.items())
    if free < 0:
        raise ValueError("pattern %s needs degree at least %d" % (pattern, d - free))

    params = ['g%d' % i for i in range(free + 1)]
    params += ['h%d_%d' % (m, i) for m, k in groups.items() for i in range(k)]
    coeffs = ['a%d' % i for i in range(d + 1)]
    R = PolynomialRing(QQ, params + coeffs)
    t = PolynomialRing(R, 't').gen()

    f = sum(R('g%d' % i) * t**i for i in range(free + 1))
    for m, k in groups.items():
        f *= (t**k + sum(R('h%d_%d' % (m, i)) * t**i for i in range(k)))**m
    f = f.list()

    graph = R.ideal([R(coeffs[i]) - f[i] for i in range(d + 1)])
    locus = graph.elimination_ideal([R(p) for p in params])
    return locus.change_ring(PolynomialRing(QQ, coeffs))


def pullback(ideal, coefficients):
    """Substitute a_i -> coefficients[i] (polynomials in a ring A) in the
    generators of an ideal of coincident_root_locus. Returns a list in A."""
    A = coefficients[0].parent()
    coefficients = [A(c) for c in coefficients]
    return [g.change_ring(A.base_ring())(*coefficients) for g in ideal.gens()]


def common_roots_conditions(f, g, k):
    """Polynomials in the coefficients of binary forms f, g (given as coefficient
    lists [a0, ..., ad]) that vanish iff f and g have >= k common roots, counted
    with multiplicity (where f = 0 identically has every root).

    The condition is that the Sylvester matrix has rank <= d1 + d2 - k.
    """
    d1, d2 = len(f) - 1, len(g) - 1
    n = d1 + d2
    rows = [[0] * i + list(f) + [0] * (d2 - 1 - i) for i in range(d2)]
    rows += [[0] * i + list(g) + [0] * (d1 - 1 - i) for i in range(d1)]
    return matrix(rows).minors(n - k + 1)
