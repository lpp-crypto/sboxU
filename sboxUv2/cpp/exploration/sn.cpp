#include "sn.hpp"


// !SECTION! Toolbox

 /**
  * Samples a random 64bits integer using Mersenne Twister from std.
  * @return The random uin64_t integer random_bits.
 */
uint64_t switching_neighbors_rand_int_64_cpp(){

    // Code for Randomness
    std::random_device rd;  /* Obtain a random seed from the OS*/
    std::mt19937_64 gen(rd());  /* Seeding the 64-bit Mersenne Twister  */
    std::uniform_int_distribution<uint64_t> dis(
        std::numeric_limits<uint64_t>::min(), 
        std::numeric_limits<uint64_t>::max()
    );

     uint64_t random_bits = dis(gen); 
    return(random_bits);
}

 /**
  * Samples n integers of size bits by sampling the minimum number of 64bits integers.
  * @param size an integer representing the size of the sampled elements.
  * @param n an integer representing the number of elements to be sampled.
  * @return The vector res of n random integers of size bits.
 */
std::vector<uint64_t> batch_rand_int_cpp(uint64_t size, uint64_t n)
{
    std::vector<uint64_t> res(n);
    uint64_t mask = 0;
    for(int i=0;i<size;i++){
        mask = mask^(1UL<<i);
    }
    uint64_t rand_int = switching_neighbors_rand_int_64_cpp();
    // Tracks the number of bits used on a random integer
    uint64_t used_bits = 0;
    // Tracks the chunk of bits we are working on
    uint64_t current_chunk = 0;
    for(int i=0;i<n;i++){
        // If we do not have enough bits on the integer, we sample a new value
        if(used_bits+size > 64){
            current_chunk = 0;
            used_bits = 0 ;
            rand_int = switching_neighbors_rand_int_64_cpp(); 
        }  
        // We take size bits located at the current chunk
        res[i] = mask & (rand_int >> (size*current_chunk));

        used_bits += size;
        current_chunk += 1;
    }
    return(res);
}

 /**
  * Returns the list of indexes of non-zero values 
  * @param s A vector of Binword
  * @return The vector of BinWord representing the indexes of non-zero values
 */

std::vector<BinWord> cpp_to_lut_coordinate(std::vector<BinWord> s){

    cpp_S_box sb = cpp_S_box(s);
    std::vector<BinWord> result;
    for(BinWord x =0; x < (1<<sb.get_input_length()); x++ ){
        if(sb[x] != 0){
            result.push_back(x);
        }
    }
    return(result);
 }  

/**
* Returns the value of the differential in a of f in x
* @param f A vector of Binword
* @param a A Binword
* @param x A Binword
* @return Returns the value of the differential in a of f in x
*/
BinWord cpp_diff_for_sn(std::vector<BinWord> f, BinWord a, BinWord x)
{
    return(cpp_oplus(f[cpp_oplus(a,x)],f[x]));
}

/**
* Returns a vector of Binword corresponding to the Bytearray
* @param v Bytearray
* @return Returns the vector of Binword corresponding to v
*/
std::vector<BinWord> cpp_from_vector( Bytearray  v){

    std::vector<BinWord> result ;
    for(BinWord i = 0; i < v.size(); i++){
        for(BinWord j = 0; j < 8; j++){
            result.push_back(1ULL & (BinWord(v[i])>>j)) ;  
        }
    }
    return(result);
}

/**
* Returns the value of the differential in a of f in x
* @param mode 'DEFAULT' or 'CUSTOM'. 'DEFAULT' uses hardcoded parameters n_eq and n_samples, while 'CUSTOM' allows to change them. 
* @param f A cpp_S_box 
* @param n_eq A uint64_t representing the maximum number of equations to reach per u
* @param n_samples A uint64_t represneting the total maximum number of equations to reach 
* @return Returns the vector of vectors of cpp_S_box representing the non-trivial switching neighbours of f
*/
 std::vector<std::vector<cpp_S_box>> cpp_non_trivial_sn (uint64_t mode,cpp_S_box f, uint64_t n_eq, uint64_t n_samples){

    
    int n = f.get_input_length();
    std::vector<std::vector<cpp_S_box>> result((1<<n)-1);
    std::vector<cpp_F2LinearSystem> E((1<<n)-1,cpp_F2LinearSystem(1<<n));

    bool valve = 1;
    while(valve){
        // We will retry until all the neighbours are APN
        bool good_neighbors = 0;

        // Adding Constraints To Remove Trivial Neighbours 
        std::vector<BinWord> temp((1<<n),0);
        for(int u = 0; u < (1<<n)-1; u++){
            // Constant Function
            temp.assign((1<<n),1);

            E[u].remove_solution(cpp_to_lut_coordinate(temp));
            // Linear functions
            for(BinWord i = 0; i < (1<<n); i++){
                // Computing an element of the linear basis
                for(BinWord x = 0; x < (1<<n); x++){
                    temp[x] = 1ULL & (x>>i);
                }
                E[u].remove_solution(cpp_to_lut_coordinate(temp));
            }
            // Removing Coordinates
            // ! TODO !
        }

        // Building The SN System
        BinWord target_u;
        std::vector<BinWord> list_sizes((f.get_lut()).size(),0);
        std::vector<BinWord> samples;
        std::vector<BinWord> ttemp(4);
        uint64_t f_size = (f.get_lut()).size();
        int c = 0;
        int c_total = 0;
        while(c < f_size-1 and c_total< n_samples){

            // Sampling 3  n bits integers by batch
            samples  = batch_rand_int_cpp(n,3);
            c_total += 1;

            // This section could be improved using tables of the differential by fixing u and finding preimages for y,y' such that y + y' = u. However, since the add_equation part seems to take the most time we leave it out for future work.

            // Computing corresponding u
            target_u = cpp_oplus(cpp_diff_for_sn(f.get_lut(),samples[0],samples[1]),cpp_diff_for_sn(f.get_lut(),samples[0],samples[2]));
            // As long as we don't have enough equations we add a new one to the corresponding constraint_list
            if (list_sizes[target_u] < n_eq){
                // Due f being APN, u = 0 means that y=x or y=x+a, so we don't have 4 '1' on the matrix   
                if (samples[1] != samples[2] and samples[1]!=cpp_oplus(samples[2],samples[0]) and samples[0]!=0 and target_u !=0){
                  ttemp[0] = cpp_oplus(samples[1],samples[0]);
                  ttemp[1] = samples[1];
                  ttemp[2] = cpp_oplus(samples[2],samples[0]);
                  ttemp[3] = samples[2];
                  E[target_u-1].add_equation(ttemp);
                list_sizes[target_u] +=1;
                }
            }
          // If the list has enough equations, we increment the counter c and increase list_size[u] so it is never modified again
          if (list_sizes[target_u] == n_eq){
              list_sizes[target_u] +=1;
              c +=1;
          }
        }

        // Computing SN
        std::vector<Bytearray> Ker;
        std::vector<BinWord> sw; 
        std::vector<BinWord> neighbour(1<<n,0);

        for(BinWord u = 1; u < (1<<n); u++){
            Ker = E[u-1].kernel_as_bytes();
            for(BinWord i = 0; i < BinWord(Ker.size()); i++){

                sw = cpp_from_vector(Ker[i]);
                // Computing the neighbour's table
                for(int k = 0; k < (1<<n); k++){
                    neighbour[k] = cpp_oplus(BinWord(u) * BinWord(sw[k]), f.get_lut()[k]) ;

                }
                if(not(cpp_is_differential_uniformity_smaller_than(neighbour,2))){
                    good_neighbors = 1;
                 }
                result[u-1].push_back(cpp_S_box(neighbour));

            }

        
        }
        // If alll functions are APN neighbours, we exit the loop
        if(not(good_neighbors)){
            valve = 0;
        }

    }
    return(result);
 }

