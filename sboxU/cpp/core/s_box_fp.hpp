#ifndef _S_BOX_FP_
#define _S_BOX_FP_

#include "../common.hpp"
#include <assert.h>

/** @class cpp_S_box_fp

    This class stores the lookup table of a vectorial function over a
    finite field \f$\mathbb{F}_p\f$ with \f$p\f$ a prime strictly
    greater than two, and provides convenient methods to interact with
    it.

    A function mapping \f$\mathbb{F}_p^t\f$ to \f$\mathbb{F}_p^u\f$ is
    represented in a "vectorized" way: every element of the input (or
    output) space is a vector of \f$t\f$ (resp. \f$u\f$) elements of
    \f$\mathbb{F}_p\f$, stored as an `FpWord` in little-endian
    convention (the coordinate of lowest weight comes first). Such a
    vector can be mapped to the integer that it represents in base
    \f$p\f$, and back, using `vec_to_int` and `int_to_vec`. The integer
    representation is used to index the lookup table.

    Storing both representations lets us rely on `std::inner_product`
    (which the compiler vectorizes aggressively) for the base-\f$p\f$ to
    integer conversion, which is the hot path when computing tables such
    as the DDT.
 */
class cpp_S_box_fp {

    private :
        BinWord input_size ;
        BinWord output_size ;
        Integer p ;
        // Contains iterated powers of p for input space, in increasing order
        std::vector<Integer> powers_in ;
        // Contains iterated powers of p for output space, in increasing order
        std::vector<Integer> powers_out ;
        // Contains the list of all input_space, in base p decomposition. The convention used is little-endian
        std::vector<FpWord> input_space ;
        // Containes the list of all output_space, in base p decomposition. The convention used is little-endian
        std::vector<FpWord> output_space;
        // Contains the list of all outputs, in base p decomposition. The convention used is little-endian
        std::vector<FpWord> lut ;

    public :
        // CONSTRUCTORS

        /** Builds an empty S-box (all sizes set to zero and all containers empty). */
        cpp_S_box_fp() : input_size(0), output_size(0), p(), powers_in(0), powers_out(0), input_space(0), output_space(0), lut(0) {}

        /** Builds an S-box from all of its precomputed components.

            This constructor performs no computation: every field is
            moved in as provided. It is meant for the case where the
            powers and spaces have already been computed (for instance
            when deriving a new S-box from an existing one) and we want
            to avoid recomputing them.

            @param _input_size The number \f$t\f$ of input coordinates.
            @param _output_size The number \f$u\f$ of output coordinates.
            @param _p The characteristic \f$p\f$ of the field.
            @param _powers_in The iterated powers of \f$p\f$ for the input space.
            @param _powers_out The iterated powers of \f$p\f$ for the output space.
            @param _input_space The base-\f$p\f$ decomposition of every input.
            @param _output_space The base-\f$p\f$ decomposition of every output.
            @param _lut The lookup table, indexed by the integer representation of the input.
         */
        cpp_S_box_fp(BinWord _input_size,
            BinWord _output_size,
            Integer _p,
            std::vector<Integer> _powers_in,
            std::vector<Integer> _powers_out,
            std::vector<FpWord> _input_space,
            std::vector<FpWord> _output_space,
            std::vector<FpWord> _lut) :
            input_size(_input_size),
            output_size(_output_size),
            p(_p),
            powers_in(std::move(_powers_in)),
            powers_out(std::move(_powers_out)),
            input_space(std::move(_input_space)),
            output_space(std::move(_output_space)),
            lut(std::move(_lut))
            {}

        /** Copy constructor. */
        cpp_S_box_fp(const cpp_S_box_fp& s) = default;

        ~cpp_S_box_fp(){}

        /** Builds an S-box from the minimal information needed to define it.

            This is the constructor most commonly used. The input size,
            output size, iterated powers and the input/output spaces are
            all deduced from the characteristic and the lookup table.

            @param _p The characteristic \f$p\f$ of the field.
            @param _lut The lookup table, where `_lut[i]` is the output
                (as an `FpWord`) on the input whose integer representation
                is `i`.
         */
        cpp_S_box_fp(Integer _p, std::vector<FpWord> _lut) : p(_p), lut(std::move(_lut)) {
            input_size = std::ceil(std::log(lut.size())/std::log(p));
            output_size = lut[0].size();

            powers_in = iterated_powers(p,input_size);
            powers_out = iterated_powers(p,output_size);

            input_space = build_input_space(p,input_size);
            output_space = build_input_space(p,output_size);
        }

        /** Reconstructs an S-box from a byte array produced by `to_bytes`.

            The byte array is expected to be laid out as follows. The
            first byte is a format marker, always equal to 0, that
            distinguishes the \f$\mathbb{F}_p\f$ serialization from the
            \f$\mathbb{F}_2\f$ one. The header then consists of three
            8-byte little-endian fields, in order: the characteristic
            \f$p\f$, the output size, and the number \f$n\f$ of entries in
            the lookup table. It is followed by a single byte giving the
            number of bytes used to encode each \f$\mathbb{F}_p\f$ value.
            The body then holds the \f$n\f$ entries; each entry stores its
            output-size many coordinates, and each coordinate is encoded
            over that many bytes in little-endian convention. This is
            exactly the format produced by `to_bytes`, of which this
            constructor is the inverse.

            @param bytes The byte array encoding an S-box.
         */
        cpp_S_box_fp(Bytearray bytes) {
            // bytes[0] is the format marker (0 for F_p); the header starts at offset 1
            Integer _p = 0;
            for (int b = 0; b < 8; b++) _p |= ((Integer)bytes[1 + b]) << (8*b);
            BinWord _output_size = 0;
            for (int b = 0; b < 8; b++) _output_size |= ((BinWord)bytes[9 + b]) << (8*b);
            Integer n = 0;
            for (int b = 0; b < 8; b++) n |= ((Integer)bytes[17 + b]) << (8*b);
            BinWord bytes_per_value = bytes[25];
            std::vector<FpWord> _lut(n, FpWord(_output_size));
            int pos = 26;
            for (int x = 0; x < n; x++)
                for (int j = 0; j < _output_size; j++){
                    BinWord v = 0;
                    for (int b = 0; b < bytes_per_value; b++)
                        v |= ((BinWord)bytes[pos++]) << (8*b);
                    _lut[x][j] = v;
                }
            *this = (n == 0) ? cpp_S_box_fp() : cpp_S_box_fp(_p, _lut);
        }

        // Getters

        /** @return The number \f$t\f$ of input coordinates. */
        BinWord get_input_size() const {return input_size;}

        /** @return The number \f$u\f$ of output coordinates. */
        BinWord get_output_size() const {return output_size;}

        /** @return The characteristic \f$p\f$ of the underlying field. */
        Integer get_p() const {return p;}

        /** @return The iterated powers of \f$p\f$ used to index the input space. */
        const std::vector<Integer>& get_powers_in() const {return powers_in;}

        /** @return The iterated powers of \f$p\f$ used to index the output space. */
        const std::vector<Integer>& get_powers_out() const {return powers_out;}

        /** @return The base-\f$p\f$ decomposition of every element of the input space. */
        const std::vector<FpWord>& get_input_space() const {return input_space;}

        /** @return The base-\f$p\f$ decomposition of every element of the output space. */
        const std::vector<FpWord>& get_output_space() const {return output_space;}

        /** @return The lookup table, indexed by the integer representation of the input. */
        const std::vector<FpWord>& get_lut() const {return lut;}

        // Operators overloading

        // Evaluation

        /** Evaluates the S-box on a given input.

            @param input The input, as an `FpWord` (a base-\f$p\f$ vector).
            @return The output of the S-box on `input`, as an `FpWord`.
         */
        FpWord operator[](const FpWord& input) const {
            return lut[vec_to_int(input,powers_in)];
        }

        /** Pointwise addition of two S-boxes over \f$\mathbb{F}_p\f$.

            The two S-boxes must share the same characteristic, input
            size and output size. The resulting S-box maps `x` to
            `(*this)[x] + s[x]`, the addition being done coordinate-wise
            modulo \f$p\f$.

            @param s The S-box to add to this one.
            @return A new S-box equal to the pointwise sum.
            @throws std::runtime_error If the sizes or characteristics do not match.
         */
        cpp_S_box_fp operator+(const cpp_S_box_fp& s) const {
            if (s.get_input_size() != input_size)
                throw std::runtime_error("Trying to add S_boxes of different input sizes");
            else if (s.get_p() != p)
                throw std::runtime_error("Trying to add S_boxes over Fp with different characteristics");
            else if (s.get_output_size() != output_size)
                throw std::runtime_error("Trying to add S_Boxes of different output sizes");
            else {
                std::vector<FpWord> new_lut(input_space.size());
                const std::vector<FpWord>& lut1 = get_lut(); const std::vector<FpWord>& lut2 = s.get_lut();
                for (int i = 0; i < new_lut.size();i++) {
                    new_lut[i] = FpWord(output_size);
                    for (int j = 0; j < output_size; j++){
                        new_lut[i][j] = (lut1[i][j]+lut2[i][j])%p;
                    }
                }
                return cpp_S_box_fp(p,new_lut);
            }
        }

        /** Pointwise subtraction of two S-boxes over \f$\mathbb{F}_p\f$.

            The two S-boxes must share the same characteristic, input
            size and output size. The resulting S-box maps `x` to
            `(*this)[x] - s[x]`, the subtraction being done
            coordinate-wise modulo \f$p\f$. The characteristic is added
            before reducing so that the result stays non-negative
            despite the unsigned representation.

            @param s The S-box to subtract from this one.
            @return A new S-box equal to the pointwise difference.
            @throws std::runtime_error If the sizes or characteristics do not match.
         */
        cpp_S_box_fp operator-(const cpp_S_box_fp& s) const {
            if (s.get_input_size() != input_size)
                throw std::runtime_error("Trying to subtract S_boxes of different input sizes");
            else if (s.get_p() != p)
                throw std::runtime_error("Trying to subtract S_boxes over Fp with different characteristics");
            else if (s.get_output_size() != output_size)
                throw std::runtime_error("Trying to subtract S_Boxes of different output sizes");
            else {
                std::vector<FpWord> new_lut(input_space.size());
                const std::vector<FpWord>& lut1 = get_lut(); const std::vector<FpWord>& lut2 = s.get_lut();
                for (int i = 0; i < new_lut.size();i++) {
                    new_lut[i] = FpWord(output_size);
                    for (int j = 0; j < output_size; j++){
                        new_lut[i][j] = (lut1[i][j] + p - lut2[i][j])%p;
                    }
                }
                return cpp_S_box_fp(p,new_lut);
            }
        }

        /** Composition of two S-boxes over \f$\mathbb{F}_p\f$.

            The output size of `s` must match the input size of this
            S-box, and both must share the same characteristic. The
            resulting S-box maps `x` to `(*this)[s[x]]`.

            @param s The inner S-box (applied first).
            @return A new S-box equal to the composition `(*this) o s`.
            @throws std::runtime_error If the sizes or characteristics do not match.
         */
        cpp_S_box_fp operator*(const cpp_S_box_fp& s) const {
            if (s.get_output_size()!=input_size)
                throw std::runtime_error("Trying to compose S_boxes but input size and output size do not match");
            if (s.get_p() != p)
                throw std::runtime_error("Trying to compose S_boxes over Fp with different characteristics");
            else {
                std::vector<FpWord> new_lut(s.get_input_space().size());
                const std::vector<FpWord>& lut1 = get_lut(); const std::vector<FpWord>& lut2 = s.get_lut();
                for (int i = 0; i<new_lut.size();i++){
                    new_lut[i] = lut1[vec_to_int(lut2[i],get_powers_in())];
                }
                return cpp_S_box_fp(p,new_lut);
            }
        }

        /** Tests two S-boxes over \f$\mathbb{F}_p\f$ for equality.

            Two S-boxes are equal when they share the same
            characteristic, input size and output size, and have
            identical lookup tables.

            @param s The S-box to compare this one to.
            @return `true` if the two S-boxes are equal, `false` otherwise.
         */
        bool operator==(const cpp_S_box_fp& s) const {
            if (p != s.get_p() || input_size != s.get_input_size() || output_size != s.get_output_size())
                return false;
            const std::vector<FpWord>& other_lut = s.get_lut();
            if (lut.size() != other_lut.size())
                return false;
            for (int x = 0; x < lut.size(); x++)
                if (lut[x] != other_lut[x])
                    return false;
            return true;
        }

        /** Builds a human-readable representation of the lookup table.

            The lookup table is rendered as a list of lists in base
            \f$p\f$, in little-endian convention; for instance
            `[[0,1],[1,0],...]`. This matches the format accepted by the
            factory when building an S-box from a lookup table.

            @return The string representation of the lookup table, or
                `"[]"` if the table is empty.
         */
        std::string content_string_repr() const {
            if (lut.size() == 0)
                return "[]";
            std::stringstream result;
            result << "[";
            for (int x = 0; x < lut.size(); x++){
                result << "[";
                for (int j = 0; j < lut[x].size(); j++){
                    result << std::dec << lut[x][j];
                    if (j + 1 < lut[x].size())
                        result << ",";
                }
                result << "]";
                if (x + 1 < lut.size())
                    result << ",";
            }
            result << "]";
            return result.str();
        }
        //
        // TODO : add inversion test, and inverse construct
        //

        /** Tests whether the S-box is a bijection.

            A necessary condition is that the input and output sizes
            match; if they do, the lookup table is scanned and the
            function is invertible if and only if it is injective.

            @return `true` if the S-box is invertible, `false` otherwise.
         */
        bool is_invertible() const {
            if (get_input_size()!=get_output_size()) return false;
            else {
                std::vector<bool> seen(pow((float)p,(float)input_size),false);
                for (const FpWord& val : lut) {
                    Integer i = vec_to_int(val,powers_in);
                    // If already seen : no injective, hence not invertible
                    if (seen[i]) return false;
                    else {
                        seen[i] = true;
                    }
                }
            return true;
            }
        }

        /** Computes the compositional inverse of the S-box.

            @return A new S-box `T` such that `T * (*this)` and
                `(*this) * T` are both the identity.
            @throws std::runtime_error If the input and output sizes
                differ, or if the S-box is not injective.
         */
        cpp_S_box_fp get_inverse() const {
            if (get_input_size()!=get_output_size()) throw std::runtime_error(
                "S-Box is not invertible : input space and output spaces sizes do not match");
            else {
                std::vector<bool> seen(pow((float)p,(float)input_size),false);
                std::vector<FpWord> inverse_lut(pow((float)p,(float)input_size));
                for (int i = 0; i < lut.size(); i++) {
                    const FpWord& val = lut[i];
                    Integer out = vec_to_int(val,powers_out);
                    // If already seen : no injective, hence not invertible
                    if (seen[out]) throw std::runtime_error("S-Box is not invertible");
                    else {
                        seen[out] = true;
                        inverse_lut[out] = int_to_vec(i,input_space);
                    }
                }
                return cpp_S_box_fp(p,inverse_lut);
            }
        }

        /** Computes the derivative of the S-box in a given direction.

            The derivative in direction `delta` is the function
            \f$x \mapsto S(x + \delta) - S(x)\f$, where both the
            addition of `delta` and the subtraction of the outputs are
            done modulo \f$p\f$.

            @param delta The direction of derivation, as an `FpWord`.
            @return A new S-box equal to the derivative in direction `delta`.
         */
        cpp_S_box_fp derivative(const FpWord& delta) const {
            std::vector<FpWord> new_lut(input_space.size());
            for (int i = 0; i < (int)input_space.size(); i++){
                // x is the actual input at integer index i (not its image)
                const FpWord& x = input_space[i];
                // compute x + delta coordinate-wise, reducing mod p
                FpWord x_delta(input_size);
                for (int j = 0; j < input_size; j++)
                    x_delta[j] = (x[j] + delta[j]) % p;
                Integer x_delta_int = vec_to_int(x_delta, powers_in);
                const FpWord& out_x_delta = lut[x_delta_int]; // S(x + delta)
                const FpWord& out_x       = lut[i];           // S(x)
                // D_delta(S)(x) = S(x+delta) - S(x) mod p, coordinate-wise
                FpWord new_out(output_size);
                for (int j = 0; j < output_size; j++)
                    new_out[j] = (out_x_delta[j] + p - out_x[j]) % p;
                new_lut[i] = new_out;
            }
            return cpp_S_box_fp(input_size,output_size,p,powers_in,powers_out,input_space,output_space,new_lut);
        }

        /** Extracts a single output coordinate of the S-box.

            @param i The index of the coordinate, where `0` is the
                coordinate of lowest weight.
            @return A new S-box, with a single output coordinate, equal
                to the `i`-th coordinate function of this S-box.
         */
        cpp_S_box_fp coordinate(const BinWord i) const {
            std::vector<FpWord> new_lut;
            for (int j = 0; j < input_space.size(); j++){
                const FpWord& out = lut[j];
                FpWord out_i;
                out_i.push_back(out[i]);
                new_lut.push_back(out_i);
            }
            return cpp_S_box_fp(p,new_lut);
        }

        /** Computes a component of the S-box.

            The component associated to `a` is the single-coordinate
            function \f$x \mapsto a \cdot S(x) \bmod p\f$, where
            \f$a \cdot S(x)\f$ is the scalar product over
            \f$\mathbb{F}_p\f$ of `a` with the output of the S-box. The
            scalar product is reduced after every term so that it stays
            in \f$[0, p)\f$.

            @param a The linear combination of output coordinates, as an
                `FpWord` of length equal to the output size.
            @return A new S-box, with a single output coordinate, equal
                to the component associated to `a`.
            @throws std::runtime_error If `a` does not have output-size
                many coordinates.
         */
        cpp_S_box_fp component(const FpWord& a) const {
            if (a.size() != output_size)
                throw std::runtime_error("Component vector does not have output_size many coordinates");
            std::vector<FpWord> new_lut(input_space.size());
            for (int x = 0; x < input_space.size(); x++){
                Integer acc = 0;
                const FpWord& out = lut[x];
                for (int j = 0; j < output_size; j++){
                    acc = (acc + a[j]*out[j]) % p;
                }
                FpWord comp;
                comp.push_back(acc);
                new_lut[x] = comp;
            }
            return cpp_S_box_fp(p,new_lut);
        }

        /** Serializes the S-box to a byte array.

            The format is self-describing so that the
            `cpp_S_box_fp(Bytearray)` constructor can reconstruct the
            S-box without any extra information. The first byte is a
            format marker, always equal to 0, that distinguishes this
            \f$\mathbb{F}_p\f$ serialization from the \f$\mathbb{F}_2\f$
            one (whose leading byte, an output bit-length, is never 0).
            The header then stores, each over 8 little-endian bytes, the
            characteristic \f$p\f$, the output size, and the number of
            entries in the lookup table; it is followed by a single byte
            giving the number of bytes used to store each
            \f$\mathbb{F}_p\f$ value. The body then stores every value of
            every entry, coordinate by coordinate, in little-endian
            convention.

            @return The byte array encoding this S-box.
         */
        Bytearray to_bytes() const {
            // number of bytes needed to store a value in [0, p-1]
            BinWord bytes_per_value = 1;
            Integer max_val = (p > 1) ? p - 1 : 1;
            while (bytes_per_value < 8 && (max_val >> (8*bytes_per_value)) != 0)
                bytes_per_value++;
            Bytearray result;
            Integer n = lut.size();
            // format marker : 0 distinguishes the F_p serialization from the F_2 one
            result.push_back(0);
            // header : p, output_size, n (each over 8 bytes) then bytes_per_value
            for (int b = 0; b < 8; b++) result.push_back((p >> (8*b)) & 0xFF);
            for (int b = 0; b < 8; b++) result.push_back((output_size >> (8*b)) & 0xFF);
            for (int b = 0; b < 8; b++) result.push_back((n >> (8*b)) & 0xFF);
            result.push_back(bytes_per_value);
            // body
            for (int x = 0; x < n; x++)
                for (int j = 0; j < output_size; j++){
                    Integer v = lut[x][j];
                    for (int b = 0; b < bytes_per_value; b++)
                        result.push_back((v >> (8*b)) & 0xFF);
                }
            return result;
        }

        /** Computes the iterated powers \f$p^0, p^1, \dots, p^{n-1}\f$.

            @param p The base.
            @param n The number of powers to compute.
            @return A vector `res` of size `n` such that `res[i]` equals
                \f$p^i\f$, in increasing order.
         */
        static std::vector<Integer> iterated_powers(Integer p, Integer n){
            std::vector<Integer> res(n);
            res[0] = 1;
            for (int i = 1; i < n; i++){
                res[i] = res[i-1]*p;
            }
            return res;
        }


        /** Builds the list of all base-\f$p\f$ decompositions of an input space.

            The integers in the range \f$[0, p^{\text{input\_size}} - 1]\f$
            are enumerated in increasing order, and each is stored as its
            base-\f$p\f$ decomposition in little-endian convention.

            @param p The characteristic.
            @param input_size The number of coordinates.
            @return A vector `res` such that `res[i]` is the base-\f$p\f$
                decomposition of the integer `i`.
         */
        static std::vector<FpWord> build_input_space(Integer p, BinWord input_size) {
            std::vector<FpWord> res = std::vector<FpWord>((Integer)pow((float)p,(float)input_size));
            res[0] = FpWord(input_size,0);
            for (int i = 1; i < res.size(); i++){
                // Increment the previous integer vector, with carry in base p
                FpWord prev = res[i-1];
                int j = 0;
                while (prev[j]==p-1){
                    prev[j]=0;
                    j++;
                }
                prev[j] += 1;
                res[i] = prev;
            }
            return res;
        }

        /** Converts an integer to its base-\f$p\f$ representation.

            The conversion is a simple lookup in a precomputed table
            (typically the input or output space of an S-box).

            @param i The integer to convert.
            @param lookup The table mapping integers to their base-\f$p\f$
                decomposition.
            @return The base-\f$p\f$ decomposition of `i`.
         */
        static FpWord int_to_vec(Integer i,const std::vector<FpWord>& lookup) {
            return lookup[i];
        }

        /** Converts a base-\f$p\f$ representation to the integer it encodes.

            The conversion is the inner product of the vector with the
            iterated powers of \f$p\f$, which the compiler vectorizes.

            @param v The base-\f$p\f$ vector to convert.
            @param powers The iterated powers of \f$p\f$ (as returned by
                `iterated_powers`).
            @return The integer represented by `v` in base \f$p\f$.
         */
        static Integer vec_to_int(const FpWord& v, const std::vector<Integer>& powers){
            return std::inner_product(v.begin(),v.end(),powers.begin(),0);
        }

};

#endif
