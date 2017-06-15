def general_surface(d):
    """returns the defining polynomial of a random plane curve of degree d"""

    mon = monomials([x1, x2, x3], [d, d, d])
    coeff = random_vector(len(mon))
    F = (coeff * vector(mon)).homogenize(x0)

    return(F)


def dual_surface(S, R1):
    """S is a principal ideal in QQ[x0,x1,x2,x3]"""
    # R1. < x0, x1, x2, x3, u0, u1, u2, u3 > = QQ['x0,x1,x2,x3,u0,u1,u2,u3']

    X = matrix(R1, [[x0], [x1], [x2], [x3]])
    U = matrix(R1, [[u0, u1, u2, u3]])
    S = S.change_ring(R1)
    J = jacobian(S.gens(), (x0, x1, x2, x3))

    IncidenceIdeal = S + ideal(J.stack(U).minors(2))
    Dual = IncidenceIdeal.saturation(ideal(x0, x1, x2, x3))[0]
    Dual = Dual.elimination_ideal((x0, x1, x2, x3))

    return Dual


# def project(S, P):
    """ returns contour image viewed from point P """
