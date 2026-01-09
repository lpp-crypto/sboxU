#include "anf.hpp"

// !SECTION! Functions on Boolean Functions


// Version naïve pour le moment
// Il existe une version optimisée en C par ex : https://www.joux.biz/algcrypt/PROGRAMS/Walsh_9-2.html

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


/////////////////////////////////////////////////


/// @brief Computes a compact representation of a quadratic function using its monomials
/// @param f a ccp_S_box
/// @return a vector of BinWord containing the representation
std::vector<BinWord> cpp_quadratic_compact_representation( const cpp_S_box &f){

    // !! CAREFUL !! 
    // We do not check is the degree is 2 for efficiency concerns

    int n = f.get_input_length();
    int m = f.get_output_length();

    // Computing the correct number of monomials
    int number_of_monomials = m * (n*(n-1)/2) ; 
    // Computing the number of BinWord needed
    int binword_number = number_of_monomials/8;
    if (number_of_monomials%8 != 0){
        binword_number += 1;
    }
    std::vector<BinWord> compact_representation(binword_number,0);
    std::vector<BinWord> anf ;

    int counter = 0;
    int index = 0;
    for(int k = 0; k < m; k++){
        anf = cpp_anf_component(f.coordinate(k));
/*         // Constant        
        if(anf[0] == 1){
            compact_representation[index] ^= 1<<counter;
        }
        //Update
        counter +=1;
        // Changing BinWord
        if(counter == 8){
            counter = 0;
            index +=1;
        }
        // Degree 1
        for(int i = 0; i < n; ++i){
            if(anf[1<<i] == 1){
                compact_representation[index] ^= 1<<counter;
                //printf(" lin %d \n", 1<<(i));
            }
            //Update
            counter +=1;
            // Changing BinWord
            if(counter == 8){
                counter = 0;
                index +=1;
            }
        } */
        // Degree 2
        for(int i = 0; i < n; ++i){
            for(int j = i+1; j < n; ++j){
                if(anf[(1<<i)+(1<<j)]== 1){
                    compact_representation[index] ^= 1<<counter;
                    //printf("quad %d \n", (1<<i) + (1<<j));
                }
                //Update
                counter +=1;
                // Changing BinWord
                if(counter == 8){
                    counter = 0;
                    index +=1;
                }
            }
        }
    }
    return(compact_representation);
}  


/// @brief Computes a lut from the compact representation of a quadratic function 
/// @param  compact_representation std::vector<BinWord>  containing the representation
/// @param n input length
/// @param m output length
/// @return a lut std::vector<BinWord> 
std::vector<BinWord> cpp_quadratic_sbox_from_compact_representation(std::vector<BinWord> compact_representation, int64_t n, int64_t m){


    std::vector<BinWord> lut(1<<n,0);
    std::vector<std::vector<BinWord>> coordinates_anf(m);
    std::vector<std::vector<BinWord>> coordinates_lut(m);

    // Computing the correct number of monomials
    int number_of_monomials = m * (n*(n-1)/2) ; 
    int n_quad = (n*(n-1)/2) ; 
    // Computing the number of BinWord needed
    int binword_number = number_of_monomials/8;
    if (number_of_monomials%8 != 0){
        binword_number += 1;
    }
    int counter = 0;
    int index = 0;
    // Computing the expanded coordinates anf
    for(int k = 0; k < m; k++){
        coordinates_anf[k].assign(1<<n,0);
/*         // Constant        
        coordinates_anf[k][0] = 1ULL & (compact_representation[index]>>counter);
        //Update
        counter +=1;
        // Changing BinWord
        if(counter == 8){
            counter = 0;
            index +=1;
        }
        // Degree 1
        for(int i = 0; i < n; i++){
            coordinates_anf[k][1<<i] = 1ULL & (compact_representation[index]>>counter);
            //Update
            counter +=1;
            // Changing BinWord
            if(counter == 8){
                counter = 0;
                index +=1;
            }
        } */
        // Degree 2
        for(int i = 0; i < n; i++){
            for(int j = i+1; j < n; j++){
                coordinates_anf[k][(1<<i)+(1<<j)] = 1ULL & (compact_representation[index]>>counter);
                //Update
                counter +=1;
                // Changing BinWord
                if(counter == 8){
                    counter = 0;
                    index +=1;
                }
            }
        }
    }
    // Computing coordinate lut using Moebius transform
    for(int k = 0; k < m; k++){
        coordinates_lut[k] = cpp_anf_component(cpp_S_box(coordinates_anf[k]));
    }
    // Computing the entire lut
    for(int x = 0; x < (1<<n); x++){
        BinWord y = 0;
        for(int k = 0; k < m; k++){
            y ^= (coordinates_lut[k][x]&1ULL) << k; 
        }
        lut[x] = y;
        y = 0;
    }

    return(lut);
}
///////////////////////////////////////////////

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
        Integer coordinate_degree = cpp_degree_component(f.component(1<<i));
        if (coordinate_degree > result)
            result = coordinate_degree;
    }
    return result;
}
