def general_surface(d):
    """returns the defining polynomial of a random plane curve of degree d"""

    mon = ((x0 + x1 + x2 + x3) ^ d).monomials()
    coeff = random_vector(len(mon))
    F = (coeff * vector(mon))

    return(ideal(F))
