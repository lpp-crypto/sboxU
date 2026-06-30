
from sage.all import *
from sboxU import *

n = 10

if __name__ == "__main__":
    for t in range(0, 20):
        print("------------")
        l = BinLinearBasis([])
        for i in range(0, 15):
            x = randint(0, 2**n)
            print(hex(x))
            l.add_to_span(x)
            span = l.span()
            y = randint(0, 2**n)
            print("{:3d} {:3d} {} | {} | {:x} in span? {} {}".format(
                l.rank(),
                rank_of_vector_set(l),
                x in span,
                l,
                y,
                l.is_in_span(y),
                y in span
            ))
