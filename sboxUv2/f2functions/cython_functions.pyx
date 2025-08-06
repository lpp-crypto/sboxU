# -*- python -*-


# !SECTION! Wrapped C++

def oplus(uint64_t x, uint64_t y):
    return cpp_oplus(x, y)


def hamming_weight(uint64_t x):
    return cpp_hamming_weight(x)


def scal_prod(uint64_t x, uint64_t y):
    return cpp_scal_prod(x, y)


def msb(uint64_t x):
    return cpp_msb(x)


def lsb(uint64_t x):
    return cpp_lsb(x)


def linear_combination(vector[uint64_t] v, uint64_t mask):
    return cpp_linear_combination(v, mask)


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

            



    
