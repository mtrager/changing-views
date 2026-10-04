"""Run with  sage -python tests/test_plane_curves.py  (or pytest)."""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from sage.all import QQ, PolynomialRing, set_random_seed

from visual_events.plane_curves import (example_curve, random_plane_curve,
                                        summary)

R = PolynomialRing(QQ, 'x0,x1,x2')


def test_examples():
    # (real flexes, real bitangents, isolated real bitangents)
    expected = {'trott': (8, 28, 0), 'two ovals': (0, 8, 4),
                'peanut': (4, 4, 2), 'nested ovals': (0, 4, 4)}
    for name, counts in expected.items():
        s = summary(example_curve(R, name))
        assert (s['real flexes'], s['real bitangents'], s['isolated real bitangents']) == counts, name
        assert s['Plücker'] and s['Klein'], name


def test_random_quartic():
    set_random_seed(3)
    s = summary(random_plane_curve(R, 4))
    assert (s['class'], s['flexes'], s['bitangents']) == (12, 24, 28)
    assert s['Plücker'] and s['Klein']


if __name__ == '__main__':
    for name, f in list(globals().items()):
        if name.startswith('test_'):
            f()
            print(name, 'ok')
