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

void cpp_sn_add_equations(cpp_S_box f, std::vector<cpp_F2LinearSystem>& E, std::vector<BinWord> indices, uint64_t n_add_eq){

    std::vector<uint64_t> size_vector(E.size(),0);
    BinWord target_u;
    std::vector<BinWord> samples;
    int n = f.get_input_length();
    std::vector<BinWord> temp(4);

    uint64_t valve = 0;
    while(valve != indices.size())
    {
        samples  = batch_rand_int_cpp(n,3);
        // Computing corresponding u
        target_u = cpp_oplus(cpp_diff_for_sn(f.get_lut(),samples[0],samples[1]),cpp_diff_for_sn(f.get_lut(),samples[0],samples[2]));
        // As long as we don't have enough equations we add a new one to the corresponding constraint_list
        if (size_vector[target_u -1] < n_add_eq and target_u != 0){
            // Due f being APN, u = 0 means that y=x or y=x+a, so we don't have 4 '1' on the matrix   
            if (samples[1] != samples[2] and samples[1]!=cpp_oplus(samples[2],samples[0]) and samples[0]!=0 and target_u !=0){
                temp[0] = cpp_oplus(samples[1],samples[0]);
                temp[1] = samples[1];
                temp[2] = cpp_oplus(samples[2],samples[0]);
                temp[3] = samples[2];
                E[target_u-1].add_equation(temp);
                size_vector[target_u-1] +=1;
            }
        }
        // If the list has enough equations, we increment the counter c and increase size_vector[u] so it is never modified again
        if (size_vector[target_u] == n_add_eq){
            valve +=1;
        }
    }
}





/**
* Returns the value of the differential in a of f in x
* @param mode 'DEFAULT' or 'CUSTOM'. 'DEFAULT' uses hardcoded parameters n_eq and n_samples, while 'CUSTOM' allows to change them. 
* @param f A cpp_S_box 
* @param n_eq A uint64_t representing the maximum number of equations to reach per u
* @param n_samples A uint64_t represneting the total maximum number of equations to reach 
* @return Returns the vector of vectors of cpp_S_box representing the non-trivial switching neighbours of f
*/
 std::vector<std::vector<cpp_S_box>> cpp_non_trivial_sn (uint64_t mode,cpp_S_box f, uint64_t n_add_eq, uint64_t n_step){

    
    int n = f.get_input_length();
    std::vector<std::vector<cpp_S_box>> result((1<<n)-1);
    std::vector<cpp_F2LinearSystem> E((1<<n)-1, cpp_F2LinearSystem(1<<n, false)); 

    ///////////////////////////////////
    // Initializing with Constraints //
    ///////////////////////////////////

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

    ////////////////////////////
    // Building The SN System //
    ////////////////////////////

    std::vector<uint64_t> indices;
    for(BinWord i = 0; i < (1<<n)-1; i++)
    {
        indices.push_back(i);
    }
    // We add n_add_eq random equations to each system
    cpp_sn_add_equations(f, E, indices, n_add_eq);

    printf("First Equations Done\n");
    /////////////////////
    // Check and Retry //
    /////////////////////

    std::vector<cpp_S_box> sn_u; 
    std::vector<uint64_t> indices_temp;
    std::vector<Bytearray> Ker;
    std::vector<BinWord> sw; 
    BinWord u; 
    std::vector<BinWord> neighbour(1<<n,0);
    int bad_neighbors = 0;

    int c = 0;
    // We repeat until all the functions are APN, i.e there are no more indices to modify
    while (indices.size() != 0)
    {
        printf("Attempt %d \n", c);
        c += 1;

        for (auto v : indices)
        std::cout << v << " ";
        std::cout << "\n";

        // We check the indices for which Kernels contain non APN functions
        for(BinWord k = 0; k < indices.size(); k++){

            // We count the number of non APN neighbours
            bad_neighbors = 0;
            // The index u correspond to the system index, i.e. the field element for switching is u+1
            u = indices[k];
            //printf("u: %ld \n", u);
            Ker = E[u].kernel_as_bytes();
            //printf("ker dim: %ld \n", Ker.size());
            for(BinWord i = 0; i < BinWord(Ker.size()); i++){
                sw = cpp_from_vector(Ker[i]);
                // Computing the neighbour's table
                for(int x = 0; x < (1<<n); x++){
                    neighbour[x] = cpp_oplus(BinWord(u+1) * BinWord(sw[x]), f.get_lut()[x]) ;
                }
                sn_u.push_back(cpp_S_box(neighbour));
               // If a neighbor is not APN, we flag this index as a bad_neighbor
                if(not(cpp_is_differential_uniformity_smaller_than(neighbour,2))){
                bad_neighbors += 1;
                }
            }

            if (bad_neighbors)
            {
                indices_temp.push_back(u);
            }
            // Otherwise, we add the good neighbors 
            else{
                result[u] = sn_u;
                sn_u.clear();
            }
            }

        // For the bad indices, we add n_step equations
        indices.clear();
        for(BinWord i = 0; i < indices_temp.size(); i++)
        {
            indices.push_back(indices_temp[i]);
        }
        indices_temp.clear();
        printf("indices length %ld \n",indices.size());
        cpp_sn_add_equations(f, E, indices, n_step);
    }
    return(result);
 }


