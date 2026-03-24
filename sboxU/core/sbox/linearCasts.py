from sage.all import Integer as sage_Integer
from sboxU.core.f2functions import oplus, i2f_and_f2i

WILDCARD = "*"


# !SECTION! General Utils


def loop_over_structure(input_lengths, masks=None):
    if masks == None:
        masks = WILDCARD * len(input_lengths)
    elif len(masks) != len(input_lengths):
        raise Exception("if masks is specified, then it must have the same length as input_length")
    elif len(input_lengths) == 0:
        return []

    if masks[-1] == WILDCARD:
        for x in range(0, 2**input_lengths[-1]):
            empty = True
            for y in loop_over_structure(input_lengths[:-1], masks[:-1]):
                empty = False
                yield y + [x]
            if empty:
                yield [x]
    else:
        empty = True
        for y in loop_over_structure(input_lengths[:-1], masks[:-1]):
            empty = False
            yield y+ [masks[0]] 
        if empty:
            # !CHECK! shouldn't this be [-1] ? 
            yield [masks[0]]
    return


def canonical_cast(x):
    if isinstance(x, (int, sage_Integer)):
        return x
    elif hasattr(x, "parent"):
        if isinstance(x.parent(), (Field)):
            return ffe_to_int(x)
        else:
            raise Exception("No canonical cast for {} of type {}".format(x, type(x)))
    elif isinstance(x, (list, tuple)):
        return [canonical_cast(x_i) for x_i in x]
    raise Exception("No canonical cast for {} of type {}".format(x, type(x)))


# !SECTION! Cartesian Products


class CastFromF2Product:
    """Encapsulates a mapping from a cartesian product of vector spaces F_2^n_i into the modular ring of integers modulo 2^{\\sum_i n_i}, the idea being that the later is the actual input for the lookup query inside an S_box instance.

    """
    def __init__(self, input_lengths, name=None):
        self.n_inputs = len(input_lengths)
        self.input_lengths = input_lengths
        self.total_input_length = sum(input_lengths)
        # setting up name
        if name == None:
            self.name = "λ_{"
            for i in input_lengths:
                self.name += "{:d},".format(i)
            self.name = "{} -> {}".format(
                self.name[:-1],
                self.total_input_length
            )
            self.name += "}"
        else:
            self.name = name
        # setting up parsing
        self.input_masks  = [ (1<<x) - 1 for x in input_lengths ]
        self.input_shifts = [ sum(input_lengths[0:i])
                              for i in range(0, self.n_inputs) ]
        self.luts = [[0 for x in range(0, 2**self.input_lengths[i])]
                     for i in range(0, self.n_inputs)]
        for i, n in enumerate(input_lengths):
            for x in range(0, 2**n):
                self.luts[i][x] = x << self.input_shifts[i]

    def is_valid_input(self, x):
        return (len(x) == self.n_inputs)

                
    def __call__(self, x):
        if self.n_inputs == 1:
            if isinstance(x, (int, sage_Integer)):
                return self.luts[0][x]
            elif isinstance(x, (list, tuple)):
                if len(x) == 1:
                    return self.__call__(x[0])
                else:
                    raise Exception("wrong input length: got {}, expected 1".format(len(x)))
        elif isinstance(x, (list, tuple)):
            y = 0
            for i, x_i in enumerate(x):
                if (x_i & self.input_masks[i]) != x_i:
                    raise Exception("Input x_{}={} is too large!".format(i, x_i))
                else:
                    y = oplus(y, self.luts[i][x_i])
            return y
        else:
            raise Exception("wrong input: {}".format(x))
            
    

    def __str__(self):
        return self.name


    def __iter__(self):
        for x in loop_over_structure(self.input_lengths):
            yield x


    def __len__(self):
        return self.n_inputs
    
            
    def structure(self, masks):
        for x in loop_over_structure(self.input_lengths, masks=masks):
            yield tuple(x)



            
class CastToF2Product:
    """Encapsulates a mapping from an integer (as used inside an S_box lookup table)
    back into a tuple of smaller integers, each representing an element of F_2^{n_i}.

    This is the inverse of CastFromF2Product: given an integer whose bit-decomposition
    concatenates several field elements, it extracts each component.

    """
    def __init__(self, output_lengths, name=None):
        self.n_outputs = len(output_lengths)
        self.output_lengths = output_lengths
        self.total_output_length = sum(output_lengths)
        # setting up name
        if name == None:
            self.name = "λ_{"
            self.name += "{:d} -> ".format(self.total_output_length)
            for i in output_lengths:
                self.name += "{:d},".format(i)
            self.name = self.name[:-2] + "}"
        else:
            self.name = name
        # setting up inner LUTs
        self.output_masks  = [ (1<<x) - 1 for x in output_lengths ]
        self.luts = [[0 for x in range(0, 2**self.total_output_length)]
                     for i in range(0, self.n_outputs)]
        for x in range(0, 2**self.total_output_length):
            y = x
            for i in range(0, self.n_outputs):
                self.luts[i][x] = y & self.output_masks[i]
                y = y >> self.output_lengths[i]


    def __call__(self, x):
        if not isinstance(x, (int, sage_Integer)):
            raise Exception("wrong argument for F2LinearRange: {}".format(x))
        else:
            result = [0 for k in range(0, self.n_outputs)]
            for k in range(0, self.n_outputs):
                result[k] = self.luts[k][x]
            return tuple(result)

        
    def __str__(self):
        return self.name



# !SECTION! Finite Field Conversions



class CastFromF2n:
    """Encapsulates a mapping from a finite field of size 2**n into the modular ring of integers modulo 2^{\\sum_i n_i}, the idea being that the later is the actual input for the lookup query inside an S_box instance.

    """
    def __init__(self, gf, name=None):
        self.gf = gf
        i2f, f2i = i2f_and_f2i(gf)
        self.inner = f2i
        # setting up name
        if name == None:
            self.name = "λ_{" + "GF({})".format(gf.cardinality()) + "}"
        else:
            self.name = name


    def is_valid_input(self, x):
        if hasattr(x, "parent"):
            return (x.parent() == self.gf)
        else:
            return False

                
    def __call__(self, x):
        return self.inner(x)
    

    def __str__(self):
        return self.name


    def input_space(self):
        for x in range(0, self.gf.cardinality()):
            yield self.mapping(x)


    def __len__(self):
        return self.gf.cardinality()



class CastToF2n:
    """Encapsulates a mapping from a finite field of size 2**n into the modular ring of integers modulo 2^{\\sum_i n_i}, the idea being that the later is the actual input for the lookup query inside an S_box instance.

    """
    def __init__(self, gf, name=None):
        self.gf = gf
        i2f, f2i = i2f_and_f2i(gf)
        self.inner = i2f
        # setting up name
        if name == None:
            self.name = "λ^{-1}_{" + "GF({})".format(gf.cardinality()) + "}"
        else:
            self.name = name

                
    def __call__(self, x):
        return self.inner(x)
    

    def __str__(self):
        return self.name


def casts_from_field(gf):
    return CastFromF2n(gf), CastToF2n(gf)

