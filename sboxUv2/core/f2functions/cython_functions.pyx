# -*- python -*-


# !SECTION! Wrapped C++

# !SUBSECTION! Bit-fiddling

def oplus(BinWord x, BinWord y):
    return cpp_oplus(x, y)


def hamming_weight(BinWord x):
    return cpp_hamming_weight(x)


def scal_prod(BinWord x, BinWord y):
    return cpp_scal_prod(x, y)


def msb(BinWord x):
    return cpp_msb(x)


def lsb(BinWord x):
    return cpp_lsb(x)


 
def linear_combination(std_vector[BinWord] v, BinWord mask):
    return cpp_linear_combination(v, mask)


def rank_of_vector_set(std_vector[BinWord] l):
    return cpp_rank_of_vector_set(l)
    

        
# !SECTION! Convenient XOR abstractions

def xor(*args):
    result = 0
    for x in args:
        if isinstance(x, int):
            result = oplus(x, result)
        else:
            for y in x:
                if isinstance(y, int):
                    result = oplus(y, result)
                else:
                    raise Exception("Trying to XOR a strange type ({})".format(type(x)))
    return result

            



    
