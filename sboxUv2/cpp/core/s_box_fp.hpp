#ifndef _S_BOX_FP_
#define _S_BOX_FP_

#include "../sboxU.hpp"

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
        cpp_S_box_fp() : input_size(0), output_size(0), p(), powers_in(0), powers_out(0), input_space(0), lut(0) {}
        cpp_S_box_fp(BinWord _input_size, 
            BinWord _output_size,
            Integer _p, 
            std::vector<Integer> _powers_in,
            std::vector<Integer> _powers_out,
            std::vector<FpWord> _input_space, 
            std::vector<FpWord> _output_space,
            std::vector<FpWord> _lut) :
            input_size(_input_size), output_size(_output_size), p(_p), powers_in(_powers_in), powers_out(_powers_out), input_space(_input_space), output_space(_output_space), lut(_lut) 
            {}
        void destruct(){
            lut.clear();
            lut.shrink_to_fit();
        }
        // This is the minimal arguments that it can take, without that it doesn't specify a SBox over Fp
        cpp_S_box_fp(Integer _p, std::vector<FpWord> _lut) : p(_p), lut(_lut) {

            input_size = std::ceil(std::log(lut.size())/std::log(p));
            output_size = std::ceil(lut[0].size());

            powers_in = iterated_powers(p,input_size);
            powers_out = iterated_powers(p,output_size);
            
            input_space = build_input_space(p,input_size);
            output_space = build_input_space(p,output_size);
        }

        // Getters
        BinWord get_input_size() const {return input_size;}

        BinWord get_output_size() const {return output_size;}

        Integer get_p() const {return p;}

        const std::vector<Integer>& get_powers_in() const {return powers_in;}

        const std::vector<Integer>& get_powers_out() const {return powers_out;}

        const std::vector<FpWord>& get_input_space() const {return input_space;}

        const std::vector<FpWord>& get_output_space() const {return output_space;}

        const std::vector<FpWord>& get_lut() const {return lut;}

        // Operators overloading

        // Evaluation
        FpWord operator[](const FpWord& input) const {
            return lut[vec_to_int(input,powers_in)];
        }

        // Pointwise addition
        cpp_S_box_fp operator+(const cpp_S_box_fp& s) const {
            if (s.get_input_size() != input_size)
                throw std::runtime_error("Trying to add S_boxes of different input sizes");
            else if (s.get_p() != p)
                throw std::runtime_error("Trying to add S_boxes over Fp with different characteristics");
            else if (s.get_output_size() != output_size)
                throw std::runtime_error("Trying to add S_Boxes of different output sizes");
            else {
                std::vector<FpWord> new_lut((Integer)pow((float)p,(float)input_size));
                const std::vector<FpWord>& lut1 = get_lut(); const std::vector<FpWord>& lut2 = s.get_lut();
                for (int i = 0; i < lut.size();i++) {
                    new_lut[i] = FpWord(input_size);
                    for (int j = 0; j < input_size; j++){
                        new_lut[i][j] = (lut1[i][j]+lut1[i][j])%p;
                    } 
                }
                return cpp_S_box_fp(p,new_lut); 
            }  
        }

        // Composition
        cpp_S_box_fp operator*(const cpp_S_box_fp& s) const {
            if (s.get_output_size()!=input_size) 
                throw std::runtime_error("Trying to compose S_boxes but input size and output size do not match");
            if (s.get_p() != p)
                throw std::runtime_error("Trying to compose S_boxes over Fp with different characteristics");
            else {
                std::vector<FpWord> new_lut((Integer)pow((float)p,(float)s.get_input_size()));
                const std::vector<FpWord>& lut1 = get_lut(); const std::vector<FpWord>& lut2 = s.get_lut();
                for (int i = 0; i<lut2.size();i++){
                    new_lut[i] = lut1[std::inner_product(lut2[i].begin(),lut2[i].end(),s.get_powers_in().begin(),0)];
                }
                return cpp_S_box_fp(p,new_lut);
            }        
        }
        //
        // TODO : add inversion test, and inverse construct
        //

        // Checks if the S-box is invertible
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

        // Returns the derivative of the S-Box, namely x -> S(x+delta) - S(x)
        cpp_S_box_fp derivative(const FpWord& delta) const {
            Integer delta_int = vec_to_int(delta,powers_in);
            std::vector<FpWord> new_lut(pow((float)p,(float)input_size));
            for (int i = 0; i < input_space.size(); i++){
                Integer j = (j+delta_int)%p;
                const FpWord& out_i = lut[i];
                const FpWord& out_j = lut[j];
                Integer out_i_int = vec_to_int(out_i,powers_out);
                Integer out_j_int = vec_to_int(out_j,powers_out);
                Integer new_out_int = (out_j_int - out_i_int)% p;
                FpWord new_out = int_to_vec(new_out_int,output_space);
                new_lut[i] = new_out;
            }
            return cpp_S_box_fp(input_size,output_size,p,powers_in,powers_out,input_space,output_space,new_lut);
        }

        // Returns the i-th coordinate function of the SBox 
        cpp_S_box_fp coordinate(const BinWord i) const {
            std::vector<FpWord> new_lut(pow((float)p,(float)input_size));
            for (int j = 0; j < input_space.size(); j++){
                const FpWord& out = lut[j];
                FpWord new_out{lut[j][i]};
                new_lut[i] = new_out;
            }
            return cpp_S_box_fp(p,new_lut);
        }

        static std::vector<Integer> iterated_powers(Integer p, Integer n){
            std::vector<Integer> res(n);
            res[0] = 1;
            for (int i = 1; i < n; i++){
                res[i] = res[i-1]*p;
            }
            return res;
        }


        // Returns the list of base p decomposition of integers in range [0,...,p^input_size - 1]
        // The integers are stored in increasing order, in little endian convention
        static std::vector<FpWord> build_input_space(Integer p, BinWord input_size) {
            std::vector<FpWord> res = std::vector<FpWord>((Integer)pow((float)p,(float)input_size));
            res[0] = FpWord(input_size,0);
            for (int i = 1; i < res.size(); i++){
                // Increment the previous integer vector
                FpWord prev = res[i-1];
                int j = 0;
                while (prev[j]==1){
                    prev[j]=0;
                    j++;
                }
                prev[j] = 1;
                res[i] = prev;
            }
            return res;
        }
        
        // Converts integer to its base p representation, using a lookup table
        static FpWord int_to_vec(Integer i,const std::vector<FpWord>& lookup) {
            return lookup[i];
        }
        // Converts a base p representation to the equivalent integer, using an inner product
        static Integer vec_to_int(const FpWord& v, const std::vector<Integer>& powers){
            return std::inner_product(v.begin(),v.end(),powers.begin(),0);
        }

};

#endif
