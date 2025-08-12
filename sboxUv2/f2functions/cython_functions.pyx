# -*- python -*-


# !SECTION! Wrapped C++

# !SUBSECTION! Bit-fiddling

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


 
def linear_combination(cpp_vector[uint64_t] v, uint64_t mask):
    return cpp_linear_combination(v, mask)


def rank_of_vector_set(cpp_vector[uint64_t] l):
    return cpp_rank_of_vector_set(l)
    

# !SUBSECTION! The Linear_basis class


cdef class Linear_basis:
    def __init__(self, cpp_vector[uint64_t] l):
        self.cpp_lb = new cpp_Linear_basis(l)

        
    def __iter__(self):
        for x in self.cpp_lb[0].get_basis():
            yield x

            
    def __len__(self):
        return self.cpp_lb[0].rank()


    def __str__(self):
        result = "("
        for x in self.cpp_lb[0].get_basis():
            result += "{:x}, ".format(x)
        return result[:-2] + ")"
        
    
    def rank(self):
        return self.cpp_lb[0].rank()

    
    def add_to_span(Linear_basis self, uint64_t x):
        self.cpp_lb[0].add_to_span(x)
        

    def span(self):
        return self.cpp_lb[0].span()
    

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

            



    
