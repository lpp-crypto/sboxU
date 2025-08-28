#ifndef _S_BOX_FP_
#define _S_BOX_FP_

#include "../sboxU.hpp"
#include <cmath>

class cpp_S_box_fp {

    private :
        BinWord input_size ;
        BinWord output_size ;
        Integer p ;
        // Contains iterated powers of p
        std::vector<Integer> powers ;
        // Contains the list of all inputs, in base p decomposition. The convention used is little-endian
        std::vector<std::vector<Integer>> inputs ;
        // Contains the list of all outputs, in base p decomposition. The convention used is little-endian
        std::vector<std::vector<Integer>> values ;

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
        cpp_S_box_fp() : input_size(0), output_size(0), p(), powers(0), inputs(0), values(0) {}
        cpp_S_box_fp(BinWord _input_size, BinWord _output_size, Integer _p, std::vector<Integer> _powers, std::vector<std::vector<Integer>> _inputs, std::vector<std::vector<Integer>> _values) :
            input_size(_input_size), output_size(_output_size), p(_p), powers(_powers), inputs(_inputs), values(_values) 
            {}
        // This is the minimal arguments that it can take, without that it doesn't specify a SBox over Fp
        cpp_S_box_fp(Integer _p, std::vector<std::vector<Integer>> _values) : p(_p), values(_values) {
            input_size = std::ceil(std::log(values.size())/std::log(p));
            output_size = std::ceil(std::log(values[0].size())/std::log(p));

            powers = std::vector<Integer>(input_size,1);
            for (int i = 1 ; i < powers.size(); i++) {
                powers[i] = powers[i-1]*p;
            }
            
            inputs = enumerate_input_space(p,input_size);
        }

        // Returns the list of base p decomposition of integers in range [0,...,p^input_size - 1]
        // The integers are stored in increasing order, in
        static std::vector<std::vector<Integer>> enumerate_input_space(Integer p, BinWord input_size) {
            std::vector<std::vector<Integer>> res = std::vector<std::vector<Integer>>((Integer)pow((float)p,(float)input_size));
            res[0] = std::vector<Integer>(input_size,0);
            for (int i = 1; i < res.size(); i++){
                // Increment the previous integer vector
                std::vector<Integer> prev = res[i-1];
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