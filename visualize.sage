def plot_curve(F, type=""):
    """ F is a homogenous ideal in QQ[x0,x1,x2] """
    v = F.gens()[0].variables()
    F = F.subs({v[0]: 1}).change_ring(QQ[v[1], v[2]])

    if (type=="surf"):
        s = singular.eval
        s('LIB "surf.lib";')
        s("ring r = 0,(" + str(v[1]) + "," + str(v[2]) + "), dp; ")
        s("plot(" + str(F.gens()[0]) + ")")

    g = F.plot(linewidth=2)

    return g


def plot_lines(lines):

    g = Graphics()
    for l in lines:
        l0 = l[0]
        line_equation = ideal(l0[0] +l0[1] *x1 +l0[2] *x2).change_ring(RR[x1, x2])
        g += line_equation.plot(color='red')

    return g


def plot_points(points_dict):

    g = Graphics()
    for p in points_dict:
        g += point((p['x1'], p['x2']), color='red', size=20)

    return g


def plot_surf(F):
    """ F is a homogeneous polynomial in x0,x1,x2,x3 """
    G = F.subs(x0=1)
    s = singular.eval
    s('LIB "surf.lib";')
    s("ring r = 0,(x1,x2,x3),dp;")
    s("surfer(" + str(G) + ")")
