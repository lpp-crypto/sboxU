# -*- python -*-



from sboxUv2.core.f2functions cimport *

from sboxUv2.core.f2functions import ffe_to_int, to_bin, from_bin, i2f_and_f2i


from sage.all import Integer as sage_Integer
from sage.all import ceil, floor


# imports needed to test the input type in the Sb factory
from sage.all import Polynomial 
from sage.crypto.sboxes import SBox as sage_SBox
from sage.rings.polynomial.multi_polynomial_element import MPolynomial

# !SECTION! Helpers


# !SUBSECTION! Functions

cdef cpp_S_box pyx_add_sboxes(cpp_S_box s, cpp_S_box t):
    """A wrapper for the cpp_S_box overloaded operator +."""
    return s.add(t)

cdef cpp_S_box pyx_mul_sboxes(cpp_S_box s, cpp_S_box t):
    """A wrapper for the cpp_S_box overloaded operator *."""
    return s.mul(t)


def new_sbox_name() -> bytes:
    """Returns a unique name that can be given to an S-box.

    It uses the module variable `sboxU_SBOXES_COUNTER` to achieve this goal by incrementing it each time it is used.

    Returns:
        A bytearray corresponding to the next unique S_box name.
    """
    global sboxU_SBOXES_COUNTER
    result = bytes("S_{}".format(sboxU_SBOXES_COUNTER).encode("UTF-8"))
    sboxU_SBOXES_COUNTER += 1
    return result


    
# !SUBSECTION! Global variables 

# global variable used to give unique names to S-boxes
cdef BinWord sboxU_SBOXES_COUNTER = 0


# !SECTION! The S_box class

cdef class S_box:
    """The S_box class stores the lookup table of an vectorial boolean function, and provides useful methods to interact with it.

    Objects of this class should be initialized using the :py:func:Sb function.

    """
    # "cdef" attributes are declared in the .pxd file
    
    # !SUBSECTION! Initialization and destruction

 
    def __init__(self, name=None):
        self.rename(name)
        self.input_cast = []
        self.output_cast = []

    
    def __dealloc__(self):
        self.cpp_sb[0].destruct()
        free(self.cpp_sb)

        
    # !SUBSECTION! Dealing with the name
    
    def rename(self, name):
        if name == None:
            self.cpp_name = new_sbox_name()
        elif isinstance(name, bytes):
            self.cpp_name = name
        elif isinstance(name, str):
            self.cpp_name = name.encode("UTF-8")
        else:
            raise NotImplemented("trying to give invalid name to S_box: {}".format(name))

        
    # !SUBSECTION! Python built-in methods
    
    def __add__(S_box self, _s) -> S_box:
        """Addition in F_2 (i.e., XOR).

        Args:
            _s: the S_box to add to the current one. Must be an S_boxable type.
        
        Returns:
            An `S_box` instance whose output is the XOR of `self` and `_s`. 
        """
        s = Sb(_s)
        if len(s) != len(self):
            raise Exception("Trying to add S_boxes of different lengths:\n{}\n{}".format(self, s))
        # below, the [0] is used to follow the S_box.cpp_sb pointer
        name = self.cpp_name + b"+" + s.name()
        result = S_box(name)
        (<S_box>result).set_inner_sbox(pyx_add_sboxes(self.cpp_sb[0], (<S_box>s).cpp_sb[0]))
        return result

    
    def __call__(self, x) -> BinWord:
        """Querying the S-box on a specific input.

        Unlike __getitem__, the input does not have to be an integer; however, it needs to be a of a type that this S_box isntance can cast to an integer. The integer obtained by querying the lookup is then cast to another type using `self.output_cast`.

        Because of the logic related to casting, it is slower than __getitem__.
        
        Args:
            x: a valid input for the cast `self.input_cast`.
        
        Returns:
            The result of calling this S-box on the input of `x`, and then casting the result to the correct type.
        """
        if isinstance(x, (int, sage_Integer)):
            return self.output_cast(self.cpp_sb.brackets(<BinWord>x))
        elif self.input_cast != None:
            for c in self.input_cast:
                if c.is_valid_input(x):
                    return self.output_cast[0](self.cpp_sb.brackets(c(x)))
            raise Exception("Could not cast input of type {} to an integer using {}".format(type(x), self.input_cast))
        else:
            raise Exception("A cast able to process {} must be specified".format(type(x)))

    # !CONTINUE! put basic casts to a dedicated file, and then add the logic to add those to S_boxes


    def __eq__(self, s) -> bool:
        if len(s) != self.cpp_sb.size():
            return False
        else:
            for x in self.input_space():
                if s[x] != self[x]:
                    return False
        return True


    def __ne__(self, s) -> bool:
        return not self.__eq__(s)

        
    def __getitem__(self, BinWord x) -> BinWord:
        """Querying the S-box on a specific integer.
        
        Args:
            x: an integer whose binary representation corresponds to the input on which to query the S-box.
        
        Returns:
            The result of calling this S-box on the input of `x`.
        """
        return self.cpp_sb.brackets(x)


    def __len__(self) -> int:
        """Returns:
            The number of entries in the lookup table of this S_box.
        """
        return self.cpp_sb.size()

    
    def __str__(self) -> str:
        return self.cpp_sb.content_string_repr().decode("UTF-8")


    def __rich_str__(self) -> str:
        if self.get_input_length() == 0:
            return "[bold][[/] [red]∅[/] [bold]][/]"
        else:
            result = "([bold blue]{:2d}[/],[bold bright_black]{:2d}[/]) [bright_black]{}[/]\n".format(
                self.get_input_length(),
                self.get_output_length(),
                self.cpp_name.decode("UTF-8"),
            )
            # first row
            if self.get_input_length() <= 4:
                n_cols = 2**self.get_input_length()
                n_rows = 1
            else:
                n_cols = 16
                n_rows = 2**(self.get_input_length() - 4)
            word_size = max(int(ceil(self.get_output_length() / 4)), 2)
            word_format = " {:" + str(word_size) + "x}"
            result += "     [black]"
            for x in range(0, n_cols):
                result += " " * (word_size-2) + " •{:x}".format(x)
            result += "[/]\n"
            # actual rows
            for i in range(0, n_rows):
                result += "[{}]  {:2x}•".format(
                    "black" if i % 2 == 0 else "blue",
                    i
                )
                
                for j in range(0, n_cols):
                    result += word_format.format(
                        self.cpp_sb.brackets(i*n_cols + j)
                    )
                result += "[/]\n"
            return result
            

    
    def __iter__(self) -> BinWord:
        for x in self.input_space():
            yield self[x]

    
    def __hash__(self):
        return hash(self.to_bytes())


    def __pow__(self, d, modulo) -> S_box:
        """Composing the S_box with itself (or with its inverse).

        Args:
            d: the number of times that the function should be iterated.

        Returns:
            If `d` is equal to 0, returns the identity S_Box. If it is a positive integer, returns the S_box corresponding to d iterations of the current S_box. If it is negative, does the same but for its inverse (throws an exception if the S_box is not invertible).
        
        The `modulo` argument is needed by the python syntax for the __pow__ function. An error will be thrown if it is set.
        """
        if modulo != None:
            raise NotImplemented("why are you using a modulo (second pow argument) here?")
        elif d == 0:
            return identity_S_box(len(self))
        elif d == 1:
            return self
        elif d == -1:
            return self.inverse()
        else:
            result = S_box(name=self.name() + b"**" + str(d).encode("UTF-8"))
            (<S_box>result).set_inner_sbox(cpp_S_box(<std_vector[BinWord]> list(self.input_space())))
            if d >= 0:
                for i in range(0, d):
                    (<S_box>result).cpp_sb[0] = pyx_mul_sboxes((<S_box>self).cpp_sb[0],
                                                              (<S_box>result).cpp_sb[0])
                return result
            elif self.is_invertible():
                inv = self.inverse()
                for i in range(0, -d):
                    (<S_box>result).cpp_sb[0] = pyx_mul_sboxes((<S_box>inv).cpp_sb[0],
                                                              (<S_box>result).cpp_sb[0])
                return result
            else:
                raise Exception("Trying to compute the negative power of a non-bijective function")


    def __lmul__(self):
        return NotImplemented("Only right composition is implemented for objects of class S_box")

    
    def __mul__(S_box self, _s) -> S_box:
        """The composition operator.

        Args:
            _s: an S_boxable object with which the current S_box must be right-composed.

        Returns:
            An S_box corresponding to the function "F o _s", where F is the current S-box, and _s is the input to the function.
        """
        s = Sb(_s)
        if self.get_input_length() != s.get_output_length():
            raise Exception("Trying to compose S_boxes of incompatible lengths:\n{}\n{}".format(self, s))
        # below, the [0] is used to follow the S_box.cpp_sb pointer
        name = self.cpp_name + "◦".encode("UTF-8") + s.name()
        result = S_box(name)
        (<S_box>result).set_inner_sbox(pyx_mul_sboxes(self.cpp_sb[0], (<S_box>s).cpp_sb[0]))
        return result

    
    # !SUBSECTION! Functions dealing with input/output sizes

    
    def get_input_length(self) -> int:
        return self.cpp_sb.get_input_length()

    
    def input_space_size(self) -> int:
        return self.cpp_sb.input_space_size()


    def input_space(self):
        return range(0, self.cpp_sb.input_space_size())

    
    def get_output_length(self) -> int:
        return self.cpp_sb.get_output_length()


    def output_space_size(self) -> int:
        return self.cpp_sb.output_space_size()

    
    def output_space(self):
        return range(0, self.cpp_sb.output_space_size())

    
    # !SUBSECTION! Basic accessors
    
    def lut(self) -> list:
        return self.cpp_sb.get_lut()


    cdef set_inner_sbox(S_box self, cpp_S_box s):
        if self.cpp_sb:
            del self.cpp_sb
        self.cpp_sb = new cpp_S_box()
        self.cpp_sb[0] = s

    
    def to_bytes(self) -> bytes:
        return bytes(self.cpp_sb.to_bytes())

    
    def name(self) -> bytes:
        return self.cpp_name

    
    # !SUBSUBSECTION! Components and coordinates
    
    def coordinate(S_box self, BinWord i) -> S_box:
        """Args:
            i: the index of the coordinate, where 0 is the bit of lowest weight.
        
        Returns:
            An S_box instance mapping n bits to 1 corresponding to the i-th coordinate of S.
        
        """
        assert i < self.cpp_sb.get_output_length()
        result = S_box(name=self.cpp_name + ("_{:x}".format(i)).encode("UTF-8"))
        result.set_inner_sbox(<cpp_S_box>self.cpp_sb.coordinate(<BinWord>i))
        return result
        
    
    def component(S_box self, BinWord a) -> S_box:
        """Returns:
            An S_box instance mapping n bits to 1 corresponding to the component x \mapsto S(x) \cdot a, where \cdot is the standard scalar product.
        
        """
        result = S_box(name=self.cpp_name + ("•{:x}".format(a)).encode("UTF-8"))
        result.set_inner_sbox(<cpp_S_box>self.cpp_sb.component(<BinWord>a))
        return result
        

    # !SUBSUBSECTION! Derivatives

    def derivative(S_box self, BinWord delta) -> S_box:
        """Returns:
            An S_box of the same dimension as S corresponding to its derivative in the direction delta, i.e. x \mapsto S(x+delta)+S(x).
        
        """
        result = S_box(name=("Δ_{:x} ".format(delta)).encode("UTF-8") + self.cpp_name)
        result.set_inner_sbox(<cpp_S_box>self.cpp_sb.derivative(<BinWord>delta))
        return result
        
    
    
    # !SUBSECTION! Function composition

    def is_invertible(self) -> bool:
        """Returns:
            True if the current S_box is a bijection, False otherwise.
        """
        return self.cpp_sb.is_invertible()

    
    def inverse(self) -> S_box:
        """Returns:
            An S_box instance corresponding to the compositional inverse of the current S_box.

        If the current S_box is not invertible, will probably crash.
        """
        name = self.cpp_name + b"^-1"
        result = S_box(name=name)
        (<S_box>result).set_inner_sbox(<cpp_S_box>(self.cpp_sb.inverse()))
        return result


    


        
# !SECTION! The S_box_fp class

cdef class S_box_fp:
    
    # !SUBSECTION! Initialization

    def __init__(self,name=None):
        self.rename(name)

    # !SUBSECTION! Dealing with the name

    def rename(self,name):
        if name == None:
            self.cpp_name = new_sbox_name()
        elif isinstance(name, bytes):
            self.cpp_name = name
        elif isinstance(name, str):
            self.cpp_name = name.encode("UTF-8")
        else:
            raise NotImplemented("trying to give invalid name to S_box: {}".format(name))

    # !SUBSECTION! Python built-in methods

    def __add__(self, _s):
        """Pointwise addition in F_p (i.e., modular addition mod p).

        Args:
            _s: the S_box to add to the current one. Must be an S_boxable type.

        Returns:
            An `S_box_fp` instance whose output is the modular addition of `self` and `_s`.

        """
        s = Sb(_s)
        if len(s) != len(self):
            raise Exception("Trying to add S_boxes of different lengths:\n{}\n{}".format(self,s))
        name = self.cpp_name + b"+" + s.name()
        result = S_box_fp(name)
        (<S_box_fp>result).set_inner_sbox((<S_box_fp>self).cpp_sb[0]+(<S_box_fp>s).cpp_sb[0])
        return result


    def __mul__(self,_s):
        """Composition of S-Boxes in F_p

        Args:
            _s (_type_): the S_box to be composed to the right to the current one. Output size and input size must match.
        
        Returns:
            An `S_box_fp` instance whose output is the composition of `self`and `_s`.
        """      
        s = Sb(_s)
        name = self.cpp_name + "◦".encode("UTF-8") + s.name()
        result = S_box_fp(name)
        (<S_box_fp>result).set_inner_sbox((<S_box_fp>self).cpp_sb[0]+(<S_box_fp>s).cpp_sb[0])
        return result
      

    def __getitem__(self, FpWord x):
        """Querying the S-box on a specific input.

        Args:
            x (FpWord): a vector of integers representing where the S-box is queried.

        Returns:
            The result of calling this S-box on the input of `x`.
        """      
        return self.cpp_sb[0][x]

    def __len__(self):
        """Returns:
            The number of entries in the lookup table of this S_box.
        """        
        return self.cpp_sb.get_lut().size()

    def __str__(self):
        return f"""S-box over F{self.cpp_sb.get_p()} \n 
        Name : {self.cpp_name} \n 
        Input size : {self.cpp_sb.get_input_size()} \n 
        Output size" : {self.cpp_sb.get_output_size()}"""

    #! SUBSECTION! Getters dealing with the underlying cpp object

    def get_p(self):
        return self.cpp_sb.get_p()

    def get_input_size(self):
        return self.cpp_sb.get_input_size()

    def get_output_size(self):
        return self.cpp_sb.get_output_size()

    def input_space_size(self):
        return pow(self.get_p(),self.get_input_size())

    def output_space_size(self):
        return pow(self.get_p(),self.get_output_size())

    def input_space(self):
        return range(0,self.get_input_size())

    def output_space(self):
        return range(0,self.get_output_size())

    def get_name(self):
        return self.cpp_name

    def lut(self):
        return self.cpp_sb.get_lut()

    cdef set_inner_sbox(S_box_fp self, cpp_S_box_fp s):
        if self.cpp_sb:
            del self.cpp_sb
        self.cpp_sb = new cpp_S_box_fp()
        self.cpp_sb[0] = s

# !SUBSECTION! Functions from the SBox

    def coordinate(S_box_fp self, BinWord i):
        """Args:
            i: the index of the coordinate, where 0 is the Fp word of lowest weight.
        
        Returns:
            An S_box instance mapping n Fp words to 1 corresponding to the i-th coordinate of S.
        
        """
        assert i < self.cpp_sb.get_output_size()
        result = S_box(name=self.cpp_name + ("_{:x}".format(i)).encode("UTF-8"))
        (<S_box_fp>result).set_inner_sbox(<cpp_S_box_fp>self.cpp_sb.coordinate(<BinWord>i))
        return result

    def derivative(S_box_fp self, FpWord delta):
        """Args:
            i: the index of the coordinate, where 0 is the bit of lowest weight.
        
        Returns:
            An S_box_fp instance mapping n Fp words to 1 corresponding to the i-th coordinate of S.
        
        """
        result = S_box(name=("Δ_{:x} ".format(delta)).encode("UTF-8")+ self.cpp_name)
        (<S_box_fp>result).set_inner_sbox(<cpp_S_box_fp>self.cpp_sb.derivative(delta))
        return result

# !SECTION! Generating S-boxes

# !SUBSECTION! Main factory

def Sb(s, name=None, input_cast=[], output_cast=None) -> S_box:
    """Turns its input into an object of the S_box class.

    If it is already an S_box instance, simply returns its
    input. Otherwise, builds the lookup table, and then create the
    corresponding S_box instance.

    Args:
        s: an object of a class that can be turned into an S_box.
        name: the name to give the object. If none is provided, one will be picked using `sboxU_SBOXES_COUNTER`.
        input_casts: a list of casts that are allowed for this S_box.
        output_cast: the function to apply to the integer output when querying the LUT.

    """
    if isinstance(s, S_box):
        return <S_box>s
    else:
        result = S_box(name=name)

        if isinstance(s, (bytes, bytearray)):
            # case of a sequence of bytes
            (<S_box>result).cpp_sb = new cpp_S_box(<Bytearray>s)

        elif isinstance(s, list):
            # Case of a list of entries. Could be a list numbers or a list of polynomials
            if len(s) == 0:
                (<S_box>result).cpp_sb = new cpp_S_box(<std_vector[BinWord]>[])
            elif isinstance(s[0], (int, sage_Integer)): # case of a lookup table
                (<S_box>result).cpp_sb = new cpp_S_box(<std_vector[BinWord]>s)
            elif isinstance(s[0], (MPolynomial)): # case of a list of polynomials
                field = s[0].parent().base_ring()
                if field.characteristic() == 2:
                    degree_vars = field.degree()
                    n_vars = len(s[0].parent().gens())
                    lut = [0 for x in range(0, 2**(n_vars*degree_vars))]
                    if degree_vars ==1: # case of an ANF
                        for x in range(0, len(lut)):
                    # duplicating the code from ..anf to prevent cross dependencies
                            x_bin = to_bin(x, n_vars)
                            y = 0
                            for i in range(0, len(s)):
                                y = (<int>(s[i](x_bin)) << i) | y
                            lut[x] = y
                        (<S_box>result).cpp_sb = new cpp_S_box(<std_vector[BinWord]>lut)
                    else : # case of multivariate polynomials on field extensions
                        i2f, f2i = i2f_and_f2i(field)
                        for x in range(0, len(lut)):
                            bin_x=to_bin(x,n_vars*degree_vars)
                            pol_inputs=[i2f(from_bin(bin_x[degree_vars*i:degree_vars*(i+1)])) for i in range(n_vars)]
                            y = 0
                            for i in range(0, len(s)):
                                y = (<int>f2i(s[i](pol_inputs)) << degree_vars*i) | y
                            lut[x] = y
                        (<S_box>result).cpp_sb = new cpp_S_box(<std_vector[BinWord]>lut)    
                else:
                    raise NotImplemented("we don't yet support polynomials over fields of characteristic >2")
            else:
                msg = "can't turn list of objects of type '{}' into an S_box".format(type(s[0]))
                print(msg)
                raise NotImplemented(msg)
            
        elif isinstance(s, sage_SBox):
            # case of a Sage-style SBox instance
            (<S_box>result).cpp_sb = new cpp_S_box(<std_vector[BinWord]>list(s))

        elif isinstance(s, BinLinearMap):
            # case of a BinLinearMap
            (<S_box>result).cpp_sb = new cpp_S_box(<std_vector[BinWord]>[])
            (<S_box>result).cpp_sb[0] = (<BinLinearMap>s).cpp_blm[0].get_cpp_S_box()

        elif isinstance(s, Polynomial):
            # case of a univariate polynomial
            field = s.base_ring()
            if field.characteristic() == 2:
                n = field.degree()
                i2f, f2i = i2f_and_f2i(field)
                lut = [f2i(s(i2f(x))) 
                       for x in range(0, 2**n)]
                (<S_box>result).cpp_sb = new cpp_S_box(<std_vector[BinWord]>lut)
            else:
                raise NotImplemented("we don't yet support polynomials over fields of characteristic >2")
        else:
            try:
                result = s.get_S_box()
            except:
                msg = "can't turn object of type '{}' into an S_box".format(type(s))
                print(msg)
                raise NotImplemented(msg)
        (<S_box>result).input_cast = input_cast
        if output_cast == None:
            (<S_box>result).output_cast = [lambda x : x]
        else:
            (<S_box>result).output_cast = output_cast
        return <S_box>result


# !SUBSECTION! Other basic structures

def identity_S_box(length) -> S_box:
    """Returns an S_box instance corresponding to the identity
    function, i.e. the one mapping x to itself.

    """
    return Sb(list(range(0, length)))


cdef S_box pyx_F2_trans(BinWord k, n):
    """Wrapper for the `cpp_translation` function. """
    result = S_box(name="Add_{}".format(k))
    result.cpp_sb = new cpp_S_box()
    result.cpp_sb[0] = cpp_translation(k, n)
    return result


def F2_trans(BinWord additive_cstte, field=None, bit_length=None) -> S_box:
    """Returns an S_box containing the lookup table of a simple XOR over a given field extension of F_2.

    If additive_cstte is an integer, then either `field` or `bit_length` must be set. If it is a field element, both `field` and `bit_length` will be ignored.
    
    Args:
        additive_cstte: the constant to add. Can be a field element or an integer. If an integer, then the field used must be specified.
        field: the field in which the multiplication must be made if `additive_cstte` is an integer.
        bit_length: the bit-length to use for both the input and output if `additive_cstte` is an integer.

    Returns:
        An S_box instance
    """
    if isinstance(additive_cstte, (int, sage_Integer)):
        k = additive_cstte
        if isinstance(bit_length, (int, sage_Integer)):
            n = bit_length
        elif "degree" in dir(field): # case of a field
            n = field.degree()
        else:
            inputs = {"field": field, "bit_length": bit_length}
            raise Exception("If `additive_cstte` is an integer then either `field` or `bit_length` must be speficied, instead, got {}".format(inputs))
    else: # case where the additive constant is a finite field element
        k = ffe_to_int(additive_cstte)
        n = additive_cstte.parent().degree()
    return pyx_F2_trans(k, n)
