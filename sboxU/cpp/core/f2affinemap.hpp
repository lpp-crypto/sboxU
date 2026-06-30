#ifndef _CORE_F2AFFINEMAP_
#define _CORE_F2AFFINEMAP_

#include "../common.hpp"
#include "./f2functions.hpp"
#include "./s_box.hpp"

// !SECTION! The cpp_F2AffineMap classs itself


/** @class cpp_F2AffineMap

    Efficiently implements an affine function using the image of the canonical basis vectors under its linear part.

    Also provides useful routines e.g. to evaluate the dimension of (the linear part of) its image. It consists of a linear part and a constant. During evaluation, the linear part is computed first, and the constant is then added to the result.
  
 */
class cpp_F2AffineMap
{
    
// !SUBSECTION! Attributes

private:
    /** The bit length of the expected inputs of the linear part. */
    Integer input_length;

    /** The expected maximum output length of the longest output of the linear part. */
    Integer output_length;

public:
    /** The images of the vectors of the canonical basis of F_2^{input_length} under the linear part. */
    std::vector<BinWord> image_vectors;

    /** The constant to add to the output of the linear part.  */
    BinWord cstte;

    
// !SUBSECTION! Constructors
    
    /** Basic constructor initialisinng everything to 0. */
    inline cpp_F2AffineMap() :
        image_vectors(std::vector<BinWord>(0)),
        input_length(0),
        output_length(0),
        cstte(0)
    {};
    
    
    /** Constructor using the images of the linear part to set the corresponding attribute, and to compute both input and output length (the constant being set to 0).  */
    inline cpp_F2AffineMap(const std::vector<BinWord> & _image_vectors) :
        image_vectors(_image_vectors.cbegin(), _image_vectors.cend()),
        input_length(_image_vectors.size()),
        output_length(0),
        cstte(0)
    {
        for(auto &v : image_vectors)
        {
            BinWord m = cpp_msb(v);
            if (m > output_length)
                output_length = m;
        }
        output_length ++;       // cpp_msb starts at 0
    };

    
    /** Constructor using the images of the linear part to set the corresponding attribute (and to compute both input and output length), and which sets the constant to the given value.  */
    inline cpp_F2AffineMap(const std::vector<BinWord> & _image_vectors,
                              const BinWord _cstte) :
        cpp_F2AffineMap(_image_vectors)
    {
        cstte = _cstte;
    };

    
    /** Constructor using the images of the linear part to set the corresponding attribute, and which forces the given values for input and output_length.

        Can be useful for instance if this mapping is intended to be composed with an other with a fixed input length, but the actual length of its image vectors is not guaranteed.

        The constant is set to 0.
    */
    inline cpp_F2AffineMap(
        const std::vector<BinWord> & _image_vectors,
        const Integer _input_length,
        const Integer _output_length
        ) :
        image_vectors(_image_vectors),
        input_length(_input_length),
        output_length(_output_length),
        cstte(0)
    {};



    /** Constructor using the images of the linear part and the constant to set the corresponding attributes, and which forces the given values for input and output_length.

        Can be useful for instance if this mapping is intended to be composed with an other with a fixed input length, but the actual length of its image vectors is not guaranteed.

        The constant is set to the given value.
    */
    inline cpp_F2AffineMap(
        const std::vector<BinWord> & _image_vectors,
        const Integer _input_length,
        const Integer _output_length,
        const BinWord _cstte
        ) :
        cpp_F2AffineMap(_image_vectors, _input_length, _output_length)
    {
        cstte = _cstte;
    };

    

    /** Constructor deriving te ehe value of all the attributes from the given lookup table.

        It only relies on the values in 0 and along the canonical basis; no check is performed.
    */
    cpp_F2AffineMap(const cpp_S_box & lut) ;

    
    
// !SUBSECTION! Getters and setters


    inline bool is_linear() const
    {
        return (cstte == 0);
    };
    
    
    inline Integer get_input_length() const
    {
        return input_length;
    };
    
    
    inline Integer get_output_length() const
    {
        return output_length;
    };


    inline cpp_S_box get_cpp_S_box() const
    {
        std::vector<BinWord> lut(1 << input_length, 0);
        for(BinWord x=0; x<lut.size(); x++)
            lut[x] = cstte ^ cpp_linear_combination(image_vectors, x);
        return cpp_S_box(lut, input_length, output_length);
    };

    
    inline std::vector<BinWord> get_image_vectors() const
    {
        std::vector<BinWord> result(image_vectors);
        return result;
    };

    inline BinWord get_cstte() const
    {
        BinWord result = cstte;
        return result;
    };

    

// !SUBSECTION! Overloaded operators

    /** Evaluate the affine mapping on the given input.

        @param x A BinWord corresponding to the binary vector on which the mapping is to be evaluated.

        @return The XOR of constant and the image of x under the linear part.
     */
    inline BinWord operator()(const BinWord x) const
    {
        return cstte ^ cpp_linear_combination(image_vectors, x);
    };


    /** Composes this F2AffineMap with another one.

        @param l Another F2AffineMap

        @return A new F2AffineMap A such A(x)=self(l(x)) for all relevant x.
     */
    cpp_F2AffineMap operator*(const cpp_F2AffineMap & l) const;


    /** Sums F2AffineMap with another one.

        @param l Another F2AffineMap

        @return A new F2AffineMap A such A(x)=self(x)+l(x) for all relevant x, where the addition is done in parallel over F_2.
     */    
    cpp_F2AffineMap operator+(const cpp_F2AffineMap & l) const;


    /** Sums F2AffineMap with a constant.

        @param cst A 64-bit integer

        @return A new F2AffineMap A such A(x)=cst + self(x) for all relevant x
     */    
    cpp_F2AffineMap operator+(const BinWord cst) const;


    /** Concatenates the output of this F2AffineMap with another one.

        @param l Another F2AffineMap

        @return A new F2AffineMap A such A(x || y)=self(x) || l(y) for all relevant x and y, where || denotes concatenation, so that the output of this F2AffineMap are the bits of lowest weight in the outputs of the result.
     */
    cpp_F2AffineMap operator|(const cpp_F2AffineMap & l) const;

    
    
// !SUBSECTION! Matrix operations


    /** Tests if compositional inversion is possible.

        @return True if and only if the linear part of this mapping has full rank, i.e. if input_length == output_length == rank.
     */
    bool is_invertible() const
    {
        if (input_length != output_length)
            return false;
        else
            return (cpp_rank_of_vector_set(image_vectors) == input_length);
    }

    
    /** Compositional inversion.

        @return A new F2AffineMap A such A * self=Id, Id being the identity. If the mapping is not invertible, throws an exception.
     */
    cpp_F2AffineMap inverse() const;


    /** In the linear case, computes a new linear F2AffineMap corresponding to the F_2 matrix obtained by transposing the one corresponding to this mapping.

        @return A new F2AffineMap A with a constant equal to 0 and with a matrix representation that is the transpose of the one for this mapping.
     */
    cpp_F2AffineMap transpose() const;


    /** The rank of the linear part of this mapping.

        @return An integer at lesat equal to zero corresponding to the rank of the linear part or, equivalently, to the dimension of the image affine space.
     */
    inline BinWord rank() const
    {
        return cpp_rank_of_vector_set(image_vectors);
    };
    
};



// !SECTION! Helper fuctions
cpp_F2AffineMap operator+(const BinWord cst, const cpp_F2AffineMap & A);

cpp_F2AffineMap identity_F2AffineMap(unsigned int n);

cpp_F2AffineMap cpp_F2AffineMap_from_lut(Lut & s);

cpp_F2AffineMap cpp_F2AffineMap_from_lut(cpp_S_box & s);

cpp_F2AffineMap cpp_block_diagonal_F2AffineMap(
    const cpp_F2AffineMap &A,
    const cpp_F2AffineMap &B
    );


#endif
