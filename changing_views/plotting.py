"""Plotting real plane curves, points and lines in P^2.

Two views are available:
  * 'disk': the projective plane as the unit disk, i.e. the upper hemisphere
    x0 >= 0 of the unit sphere seen from above. The whole real curve is
    visible; antipodal points of the boundary circle are identified.
  * 'affine': the affine chart x0 = 1, restricted to a window.
"""
from math import nan, sqrt

from sage.all import Graphics, circle, fast_callable, implicit_plot, point, vector


def _evaluator(f):
    if hasattr(f, 'gens') and not hasattr(f, 'variables'):  # an ideal
        f = f.gens()[0]
    return fast_callable(f, vars=f.parent().gens(), domain=float)


def _zero_set(g, view, window, **kwds):
    """Plot {g(x0, x1, x2) = 0} for a fast callable g of three arguments."""
    kwds.setdefault('plot_points', 400)
    if view == 'disk':
        def h(a, b):
            r = 1 - a * a - b * b
            return g(sqrt(r), a, b) if r >= 0 else nan
        return implicit_plot(h, (-1, 1), (-1, 1), **kwds)
    return implicit_plot(lambda a, b: g(1, a, b), window, window, **kwds)


def _frame(view):
    if view == 'disk':
        return circle((0, 0), 1, color='gray', linestyle='--', aspect_ratio=1, axes=False)
    return Graphics()


def plot_curve(f, view='disk', window=(-3, 3), color='blue', **kwds):
    """The real points of the plane curve {f = 0} (a form or principal ideal)."""
    kwds.setdefault('linewidth', 1.5)
    return _frame(view) + _zero_set(_evaluator(f), view, window, color=color, **kwds)


def plot_lines(lines, view='disk', window=(-3, 3), color='red', **kwds):
    """Lines u0*x0 + u1*x1 + u2*x2 = 0, given as real vectors u."""
    kwds.setdefault('linewidth', 0.8)
    G = Graphics()
    for u in lines:
        u = [float(c) for c in u]
        G += _zero_set(lambda x0, x1, x2: u[0] * x0 + u[1] * x1 + u[2] * x2,
                       view, window, color=color, **kwds)
    return G


def to_plane(p, view='disk'):
    """Coordinates in the picture of a real point p of P^2."""
    p = vector([float(c) for c in p])
    if view == 'disk':
        p = p / p.norm()
        if p[0] < 0:
            p = -p
        return (p[1], p[2])
    return (p[1] / p[0], p[2] / p[0])


def plot_points(points, view='disk', color='red', size=25, **kwds):
    G = Graphics()
    for p in points:
        if view == 'affine' and abs(float(p[0])) < 1e-12:
            continue  # at infinity
        G += point(to_plane(p, view), color=color, size=size, zorder=5, **kwds)
    return G
