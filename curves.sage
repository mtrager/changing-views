def general_curve(d):
    """returns the ideal of a random plane curve of degree d"""

    mon = ((x0 + x1 + x2) ^ d).monomials()
    coeff = random_vector(len(mon))
    F = (coeff * vector(mon))

    return(ideal(F))


def dual_curve(F, R1):
    """F is an equation in QQ[x0,x1,x2]"""
    # R1. < x0, x1, x2, u0, u1, u2 > = QQ['x0,x1,x2,u0,u1,u2']

    F = F.change_ring(R1)
    J = jacobian(F.gens()[0], (x0, x1, x2))

    U = matrix(R1, [[u0, u1, u2]])
    incidence_ideal = F + ideal(J.stack(U).minors(2))
    dual_curve = incidence_ideal.saturation(ideal(x0, x1, x2))[0]
    dual_curve = dual_curve.elimination_ideal((x0, x1, x2))

    return dual_curve


# subs(dict{zip(z,vals)})


def plot_curve(F):
    """ F is a homogeneous polynomial in x0,x1,x2,x3 """
    v = F.gens()[0].variables()
    F = F.subs({v[0]: 1})
    s = singular.eval
    s('LIB "surf.lib";')
    s("ring r = 0,(" + str(v[1]) + "," + str(v[2]) + "), dp; ")
    s("plot(" + str(F.gens()[0]) + ")")
