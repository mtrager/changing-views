def discriminant_polynomial(degree):
    """returns the discriminant of a polynomial of degree d
    in the variables (coefficients) a0,...,ad"""

    v=['a%d' % i for i in xrange(degree +1)]
    v.append('t')
    R=PolynomialRing(QQ, v)
    R.inject_variables()
    poly = sum([R('a%d' %i) *t ^i for i in range(degree +1)])
    res = poly.resultant(poly.derivative(t), t)
    res_ideal = ideal(res).saturation(ideal(R('a%d' %degree)))[0]

    res = res_ideal.gens()[0]
    return res


def subresultant_1(poly0, poly1):

    d0=poly0.degree()
    d1=poly1.degree()
    D=d0 +d1
    sylv1 = (poly0.sylvester_matrix(poly1).delete_rows([d1 - 1, D -1]))
    sylv1 = sylv1.delete_columns([D -1])

    return sylv1.delete_columns([D -2]).determinant(), sylv1.delete_columns([D -3]).determinant()


def subresultant_2(poly0, poly1):

    d0=poly0.degree()
    d1=poly1.degree()
    D=d0 +d1
    sylv1 = (poly0.sylvester_matrix(poly1).delete_rows([d1 - 1, d1 -2, D -1, D -2]))
    sylv1 = sylv1.delete_columns([D -1, D -2])

    c0=sylv1.delete_columns([D -3, D -4]).determinant()
    c1=sylv1.delete_columns([D -3, D -5]).determinant()
    c2=sylv1.delete_columns([D -4, D -5]).determinant()

    return c2, c1, c0


# Old stuff

def nodal_and_cuspidal_discriminant_ideals(degree):

    v=['a%d' % i for i in xrange(degree +1)]
    v.append('t')
    R=PolynomialRing(QQ, v)
    R.inject_variables()
    poly = sum([R('a%d' %i) *t ^i for i in range(degree +1)])

    poly = poly.polynomial(t)
    polyD = poly.derivative(t)

    subres0 = poly.polynomial(t).discriminant()
    subres1 = subresultant_1(poly, polyD)
    subres2 = subresultant_2(poly, polyD)

    disc_subres2 = (subres2[1] ^2 -4 *subres2[0] *subres2[2])

    sing_disc = ideal(subres0) +ideal(subres1)

    # cuspidal
    cuspidal = sing_disc +ideal(disc_subres2)
    cuspidal = cuspidal.saturation(subres2[0])[0]

    # nodal
    nodal = sing_disc.saturation(disc_subres2)[0]

    return nodal, cuspidal


def nodal_and_cuspidal_discriminant_ideals_old(degree):

    disc=discriminant_polynomial(degree)
    R=parent(disc)
    j=disc.jacobian_ideal()
    dec=j.primary_decomposition()
    nodal=[]
    cuspidal=[]
    for i in dec:
        if (i.is_prime()):
            nodal.append(i)
        else:
            cuspidal.append(i.radical())
    if len(nodal)==1:
        nodal=nodal[0]
    if len(cuspidal)==1:
        cuspidal=cuspidal[0]
    return nodal, cuspidal


def cuspidal_discriminant_ideal(degree):
    """coefficients of polynomials with one triple root"""

    disc=discriminant_polynomial(degree)
    R=parent(disc)
    j=disc.jacobian_ideal()
    dec=j.primary_decomposition()
    cuspidal=[i.radical() for i in dec if not i.is_prime()]
    if len(cuspidal)==1:
        cuspidal=cuspidal[0]
    return cuspidal


def nodal_discriminant_ideal(degree):
    """coefficients of polynomials with two double roots"""

    disc=discriminant_polynomial(degree)
    R=parent(disc)
    j=disc.jacobian_ideal()
    dec=j.primary_decomposition()
    nodal=[i for i in dec if i.is_prime()]
    if len(cuspidal)==1:
        nodal=nodal[0]
    return nodal
