def plot_surf(F):
    """ F is a homogeneous polynomial in x0,x1,x2,x3 """
    G = F.subs(x0=1)
    s = singular.eval
    s('LIB "surf.lib";')
    s("ring r = 0,(x1,x2,x3),dp;")
    s("surfer(" + str(G) + ")")
