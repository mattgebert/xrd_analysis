import math

class spacing:
    """Spacing is a class that offers crystallographic space functions.
    Specifically, for a set of lattice lengths and a set of miller indexes (i,j,k), 
    the distance between crystallographic planes is returned.
    """
    
    def hex(lattice, miller):
        a, b, c = lattice
        h, k, l = miller
        d1 = 4/3 * (h**2 + h*k + k**2)/a**2
        d2 = l**2 / c**2
        return math.sqrt(1/(d1 + d2))


    def cubic(lattice, miller):
        a, b, c = lattice
        h, k, l = miller
        d1 = h**2 + k**2 + l**2
        return a * math.sqrt(1/(d1))
