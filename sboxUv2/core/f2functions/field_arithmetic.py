from sage.all import version, Integer
import sage


# a key issue with finite field arithmetic is that functions names
# change depending on the version of SAGE. The main purpose of this
# file is to ensure that we can write version-independent functions
# interacting with finite fields.


SAGE_VERSION = tuple([int(x) for x in sage.version.version.split(".")])

if SAGE_VERSION < (9, 8):
    def i2f_and_f2i(gf):
        """Returns the functions mapping field elements to integers
        (f2i) and integers to field elements (i2f) as a pair.


        Args:
            - gf: the finite field with which we want to interact.
              Could have been obtained using e.g. GF(q)

        """
        if gf.characteristic() == 2:
            return gf.fetch_int, lambda x : x.integer_representation() 
        else:
            return gf.__call__, lambda x : Integer(x)

        
    def ffe_from_int(x, gf):
        if gf.characteristic() > 2:
            return gf(x)
        else:
            return gf.fetch_int(x)

        
    def ffe_to_int(x):
        if x == 0:         # necessary because of inconsistent casting
            return x
        elif x == 1:       # same...
            return x
        elif x.base_ring().characteristic() > 2:
            return Integer(x)
        else:
            return x.integer_representation()
        
else:
    def i2f_and_f2i(gf):
        """A Helper function to deal with finite field elements and their integer representations.

        Returns:
            A pair of functions, namely the functions mapping field elements to integers (f2i) and the one mapping integers to field elements (i2f).


        Args:
            gf: the finite field with which we want to interact. Could have been obtained using e.g. GF(q)

        """
        return gf.from_integer, lambda x : x.to_integer()


    def ffe_from_int(x, gf):
        return gf.from_integer(x)

    
    def ffe_to_int(x):
        return x.to_integer()

