#include "anf.hpp"

// !SECTION! Functions on Boolean Functions


// Naive version for now
// See  https://www.joux.biz/algcrypt/PROGRAMS/Walsh_9-2.html for an optimized version

std::vector<BinWord> cpp_anf_component( const cpp_S_box &f)
{   if (f.get_output_length()!=1){
        throw std::runtime_error("This function is for boolean functions only");
    }
    else{
        int n= f.get_input_length();
        int N = 1<<n;
        std::vector<BinWord> v = f.get_lut();

        for (int bit = 0; bit < n; ++bit) {
            for (int mask = 0; mask < N; ++mask) {
                if (mask & (1 << bit)) {
                    v[mask] ^= v[mask ^ (1 << bit)];
                }
            }
        }
        return v;
    }
}

// Returns the degree of a component
Integer cpp_degree_component(const cpp_S_box &f)
{
    if (f.get_output_length()!=1){
        throw std::runtime_error("This function is for boolean functions only");
    }
    else{
        int n= f.get_input_length();
        int N = 1<<n;
        Integer res=-1;
        std::vector<BinWord> anf= cpp_anf_component(f);
        for (int x = 0; x < N; ++x) {
            if (anf[x] == 1) {
                int w = cpp_hamming_weight(x);
                if (w > res) {
                    res = w;
                }
            }
        }
        return res;
    }
}

// Returns the number of monomials of each degree for a component
cpp_Spectrum cpp_monomial_degree_spectrum_component(const cpp_S_box &f)
{
    if (f.get_output_length()!=1)
    {
        throw std::runtime_error("This function is for boolean functions only");
    }
    else
    {
        int n= f.get_input_length();
        int N = 1<<n;
        cpp_Spectrum mon_deg_spec=cpp_Spectrum();
        std::vector<BinWord> anf= cpp_anf_component(f);
   
        for (int x = 0; x < N; ++x)
            if (anf[x]==1)
                mon_deg_spec.incr(cpp_hamming_weight(x));
        return mon_deg_spec;
    }
}

// !SECTION! Functions for Vectorial Boolean Functions

cpp_Spectrum cpp_degree_spectrum(const cpp_S_box &f)
{   
    int m= f.get_output_length();
    int M = 1<<m;
    cpp_Spectrum deg_spec=cpp_Spectrum();
    for (int u=1; u<M; ++u)
    {
        // !TODO! use multi-threading in cpp_degree_spectrum 
        deg_spec.incr(cpp_degree_component(f.component(u)));
    }
    return deg_spec;
}


Integer cpp_algebraic_degree(const cpp_S_box &f)
{
    Integer result = 0;
    for (unsigned int i=0; i<f.get_output_length(); i++)
    {
        Integer coordinate_degree = cpp_degree_component(f.coordinate(i));
        if (coordinate_degree > result)
            result = coordinate_degree;
    }
    return result;
}

bool cpp_is_degree_bigger_than_component(const cpp_S_box &f, Integer d ){
    if (f.get_output_length()!=1)
    {
        throw std::runtime_error("This function is for boolean functions only");
    }
    else
    {
        int n= f.get_input_length();
        int N = 1<<n;
        cpp_Spectrum mon_deg_spec=cpp_Spectrum();
        std::vector<BinWord> anf= cpp_anf_component(f);
   
        for (int x = 0; x < N; ++x)
            if (anf[x]==1 & cpp_hamming_weight(x) > d){
                return true;
            }
        return false;
    }
}

bool cpp_is_degree_bigger_than(const cpp_S_box &f, Integer d ){
    for (unsigned int i=0; i<f.get_output_length(); i++){
        if (cpp_is_degree_bigger_than_component(f.coordinate(i),d)){
            return true;
        }
    }
    return false;
}