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


def subresmat(p, q, i):
    mat = p.sylvester_matrix(q)
    dq = q.degree()
    mat=mat.delete_rows(list(range(dq, dq +i)))
    mat=mat.delete_rows(list(range(i)))
    mat=mat.delete_columns(list(range(i)))
    return mat


def subrescoefs(p, q, i):
    # p, q univariate
    mat = subresmat(p, q, i)
    k = mat.nrows()
    cols = range(k -1)
    return [mat.matrix_from_columns(cols +[j]).det() for j in range(k -1, k +i)]

def subdisc_coefs(p, i):
    # p, q univariate
    q = p.derivative()
    disc = subrescoefs(p,q,i)
    return []

    # Old stuff


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

## COMPUTATION OF DISCRIMINANTS

degree = 6
#order = '3+2'

v=['a%d' % i for i in xrange(degree +1)]
v.append('t')
R=PolynomialRing(QQ, v)
R.inject_variables()
poly = sum([R('a%d' %i) *t ^i for i in range(degree +1)]).polynomial(t)

poly_d = poly.derivative()

S0 = subrescoefs(poly, poly_d, 0)
S1 = subrescoefs(poly, poly_d, 1)
S2 = subrescoefs(poly, poly_d, 2)
S3 = subrescoefs(poly, poly_d, 3)


# disc

dS = load('discriminants/disc2_2')
sub_dict = dict(zip((dS.ring()).gens(),S2))
S2disc_2 = dS.substitute(sub_dict)

dS = load('discriminants/disc3_3')
sub_dict = dict(zip((dS.ring()).gens(),S3))
S3disc_3 = dS.substitute(sub_dict)

dS = load('discriminants/disc3_2')
sub_dict = dict(zip((dS.ring()).gens(),S3))
S3disc_2 = dS.substitute(sub_dict)


if order == '2':
    disc = ideal(S0).saturation(ideal(S0[0].variables()[-1]))[0]

if order == '3':
    disc = ideal(S0)+ideal(S1)+S2disc_2
    disc = disc.saturation(S2[0])[0]
    disc = disc.saturation(ideal(S0[0].variables()[-1]))[0]

if order == '2+2':
    disc = ideal(S0[0])+ideal(S1)
    disc = disc.saturation(ideal(S0[0].variables()[-1]))[0]
    disc = disc.saturation(S2disc_2)[0]
    
if order == '4':
    disc = ideal(S0)+ideal(S1)+ideal(S2)+S3disc_3
    disc = disc.saturation(ideal(S0[0].variables()[-1]))[0]
    disc = disc.saturation(S3[0])[0]
    
if order == '3+2':
    disc = ideal(S0)+ideal(S1)+ideal(S2)+S3disc_2
    disc = disc.saturation(ideal(S0[0].variables()[-1]))[0]
    disc = disc.saturation(S3disc_3)[0]
    disc = disc.saturation(S3[0])[0]

if order == '2+2+2':
    disc = ideal(S0)+ideal(S1)+ideal(S2)
    disc = disc.saturation(ideal(S0[0].variables()[-1]))[0]
    disc = disc.saturation(S3disc_2)[0]
