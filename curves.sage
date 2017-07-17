def general_curve(d):
    """returns the ideal of a random plane curve of degree d"""

    mon = ((x0 + x1 + x2) ^ d).monomials()
    coeff = random_vector(len(mon))
    F = (coeff * vector(mon))

    return(ideal(F))


def dual_curve(F):
    """F is an equation in QQ[x0,x1,x2]"""

    R2.<x0, x1, x2, u0, u1, u2, t> = QQ[]

    # points on u0x0+u1x1+u2x2=0
    P = vector([-u1, u0, 0])
    Q = vector([0, -u2, u1])

    F = F.change_ring(R2)
    poly = (F.subs(x0=((P + t * Q)[0]), x1=((P + t * Q)[1]), x2=((P + t * Q)[2]))).gens()[0]
    dual_curve = poly.resultant(poly.derivative(t), t)

    # remove extra components that come from content respect to u0
    dual_curve = R2(dual_curve / reduce(gcd, dual_curve.polynomial(u0).coefficients()))

    return dual_curve


def singular_points(F):
    """F is a homogeneous ideal in QQ[x0,x1,x2]"""

    J = jacobian(F.gens()[0], F.gens()[0].variables())

    return F + ideal(J.minors(1))


def hessian_determinant_curve(F):

    return ideal(jacobian(jacobian(F.gens()[0], (x0, x1, x2)), (x0, x1, x2)).determinant())


def tangents_at_points(F, point_dict):

    lines = []
    J = (jacobian(F.gens()[0], (x0, x1, x2))).subs(x0=1)
    for p in point_dict:
        lines.append(J.subs(x1=p['x1'], x2=p['x2']))

    return lines
