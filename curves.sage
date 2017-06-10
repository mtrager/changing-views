def general_curve(d):
    """returns the ideal of a random plane curve of degree d"""

    mon = ((x0 + x1 + x2) ^ d).monomials()
    coeff = random_vector(len(mon))
    F = (coeff * vector(mon))

    return(ideal(F))


def dual_curve(F):
    """F is an equation in QQ[x0,x1,x2]"""

    R2.< x0, x1, x2, u0, u1, u2, t > = QQ[]

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


def hessian_determinant(F):

    return ideal(jacobian(jacobian(F.gens()[0], (x0, x1, x2)), (x0, x1, x2)).determinant())


def plot_curve(F, type=""):
    """ F is a homogenous ideal in QQ[x0,x1,x2] """
    v = F.gens()[0].variables()
    F = F.subs({v[0]: 1}).change_ring(QQ[v[1], v[2]])

    if (type=="surf"):
        s = singular.eval
        s('LIB "surf.lib";')
        s("ring r = 0,(" + str(v[1]) + "," + str(v[2]) + "), dp; ")
        s("plot(" + str(F.gens()[0]) + ")")

    P = F.plot()

    return P


def find_flex_lines(F):

    flex_points = F + hessian_determinant(F)
    zero = (flexes.subs(x0=1).change_ring(QQ[x1, x2])).variety(RR)
    l=[]
    for p in zero:
        return True


### OLD STUFF ###

def dual_curve_old(F, R1):
    """F is an equation in QQ[x0,x1,x2]"""
    # R1. < x0, x1, x2, u0, u1, u2 > = QQ['x0,x1,x2,u0,u1,u2']

    F = F.change_ring(R1)
    J = jacobian(F.gens()[0], (x0, x1, x2))

    U = matrix(R1, [[u0, u1, u2]])
    incidence_ideal = F + ideal(J.stack(U).minors(2))
    dual_curve = incidence_ideal.saturation(ideal(x0, x1, x2))[0]
    dual_curve = dual_curve.elimination_ideal((x0, x1, x2))

    return dual_curve
