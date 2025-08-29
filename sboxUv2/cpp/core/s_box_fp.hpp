#ifndef _S_BOX_FP_
#define _S_BOX_FP_

#include "../sboxU.hpp"
#include <cmath>
#include <numeric>

class cpp_S_box_fp {

    private :
        BinWord input_size ;
        BinWord output_size ;
        Integer p ;
        // Contains iterated powers of p
        FpWord powers ;
        // Contains the list of all input_space, in base p decomposition. The convention used is little-endian
        std::vector<FpWord> input_space ;
        // Contains the list of all outputs, in base p decomposition. The convention used is little-endian
        std::vector<FpWord> lut ;

        // Performs fast exponentiation x^k mod p
        Integer mod_pow(Integer x, Integer k, Integer p) {
            Integer res = x;
            while (k > 0) {
                if (k % 2 == 0) {res = (res*res)%p; k=k/2;}
                else {res = (res*res*x)%p; k=k/2;}
            }
            return res;
        }
    
    public :
        // CONSTRUCTORS
        cpp_S_box_fp() : input_size(0), output_size(0), p(), powers(0), input_space(0), lut(0) {}
        cpp_S_box_fp(BinWord _input_size, 
            BinWord _output_size,
            Integer _p, 
            FpWord _powers, 
            std::vector<FpWord> _input_space, 
            std::vector<FpWord> _lut) :
            input_size(_input_size), output_size(_output_size), p(_p), powers(_powers), input_space(_input_space), lut(_lut) 
            {}
        // This is the minimal arguments that it can take, without that it doesn't specify a SBox over Fp
        cpp_S_box_fp(Integer _p, std::vector<FpWord> _lut) : p(_p), lut(_lut) {

            input_size = std::ceil(std::log(lut.size())/std::log(p));
            output_size = std::ceil(std::log(lut[0].size())/std::log(p));

            powers = FpWord(input_size,1);
            for (int i = 1 ; i < powers.size(); i++) {
                powers[i] = powers[i-1]*p;
            }
            
            input_space = build_input_space(p,input_size);
        }

        // Getters
        BinWord get_input_size() const {return input_size;}

        BinWord get_output_size() const {return output_size;}

        Integer get_p() const {return p;}

        const FpWord& get_powers() const {return powers;}

        const std::vector<FpWord>& get_input_space() const {return input_space;}

        const std::vector<FpWord>& get_lut() const {return lut;}

        // Operators overloading

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
                    new_lut[i] = lut1[std::inner_product(lut2[i].begin(),lut2[i].end(),s.get_powers().begin(),0)];
                }
                return cpp_S_box_fp(p,new_lut);
            }        
        }
        //
        // TODO : add inversion test, and inverse construct
        //

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

};

#endif