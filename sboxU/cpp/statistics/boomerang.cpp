#include "boomerang.hpp"
#include <array>



// !SECTION! BCT 
// !SUBSECTION! The BCT itself 


std::vector<Integer> cpp_bct_row(
    const cpp_S_box & s,
    const BinWord a)
{
    if (a == 0)
        return std::vector<Integer>(s.output_space_size(),s.input_space_size());

    std::vector<Integer> result(s.output_space_size(), 0);

    // Flat linked-list structure replacing vector<vector<BinWord>>.
    // For each output difference z, we maintain a singly-linked list of all
    // y1 values seen so far with that z. This avoids N heap allocations and
    // keeps all data cache-local in three contiguous arrays.
    //
    // pool[i]  : the i-th stored y1 value
    // head[z]  : index in pool of the most recently inserted y1 for bucket z
    //            (-1 if the bucket is empty)
    // next[i]  : index in pool of the previous y1 in the same bucket as pool[i]
    //            (-1 if pool[i] is the oldest entry in its bucket)
    //
    // At most N/2 entries are inserted in total (one per coset representative).
    std::vector<BinWord> pool(s.input_space_size() / 2);
    std::vector<int> head(s.output_space_size(), -1);
    std::vector<int> next(s.input_space_size() / 2, -1);
    int pool_idx = 0;

    result[0] = s.input_space_size();

    FOR_ENUMERATE_DIFFERENCE_COSETS(x,a,s.input_space_size())
    {
        BinWord y1 = s[x];
        BinWord z = y1 ^ s[x^a];

        // For each previously seen y2 with the same output difference z,
        // both (y1^y2) and (y1^y2^z) are valid output differences for the BCT.
        for(int i = head[z]; i != -1; i = next[i])
        {
            result[ y1 ^ pool[i] ] += 4;
            result[ y1 ^ pool[i] ^ z] += 4;
        }

        // Insert y1 at the head of the linked list for bucket z.
        pool[pool_idx] = y1;
        next[pool_idx] = head[z];
        head[z] = pool_idx++;

        result[z] += 2;
    }
    return result;
}

std::vector< std::vector<Integer>> cpp_bct(const cpp_S_box & s)
{
	std::vector< std::vector<Integer>> table_bct ;
        table_bct.reserve(s.input_space_size());
	for(unsigned int a = 0; a < s.input_space_size(); a++)
		table_bct.push_back(cpp_bct_row(s,a)) ;
	return table_bct ;
}


// !SUBSECTION! BCT Spectrum 


cpp_Spectrum cpp_boomerang_spectrum(const cpp_S_box & s, const unsigned int n_threads)
{
    cpp_Spectrum count;
    int threads = threads_from_size(s.input_space_size());
#pragma omp parallel for reduction(aggregateSpectrum:count) num_threads(threads)
    for( unsigned int a = 1; a < s.input_space_size(); a++)
        count.incr_by_counting(cpp_bct_row(s,a));
    count.incr_by_amount(s.input_space_size(),-(s.input_space_size()-1));
    if(count[s.input_space_size()] == 0)
        count.erase(s.input_space_size());
    return count;
}

// !SECTION! FBCT

std::vector<Integer> cpp_fbct_row(
    const cpp_S_box & s,
    const BinWord a)
{
    if (a == 0)
        return std::vector<Integer>(s.input_space_size(),s.input_space_size());

    std::vector<Integer> result(s.input_space_size(), 0);
    std::vector<std::vector<BinWord> > xor_list(s.input_space_size(), std::vector<BinWord>(0));
    // result[0] = s.input_space_size();
    // result[a] = s.input_space_size();

    FOR_ENUMERATE_DIFFERENCE_COSETS(x,a,s.input_space_size())
    {
        BinWord z = s[x] ^ s[x^a];

        for(auto x2 : xor_list[z]){
            result[ x ^ x2 ] += 4;
            result[ x ^ x2 ^ a ] += 4;
        }
        xor_list[z].push_back(x);
    }
    return result;
}

cpp_Spectrum cpp_fbct_spectrum(const cpp_S_box & s, const unsigned int n_threads)
{
    cpp_Spectrum count;
    int threads = threads_from_size(s.input_space_size());
#pragma omp parallel for reduction(aggregateSpectrum:count) num_threads(threads)
    for( unsigned int a = 1; a < s.input_space_size(); a++)
        count.incr_by_counting(cpp_fbct_row(s,a));
    count.incr_by_amount(0,-s.input_space_size()*2+2);
    return count;
}


std::vector< std::vector<Integer>> cpp_fbct(
    const cpp_S_box & s)
{
    std::vector<std::vector<Integer>> result(s.input_space_size(), std::vector<Integer>(s.input_space_size(), 0));
    result[0] = std::vector<Integer>(s.input_space_size(), s.input_space_size());
    result[1][0] = s.input_space_size();
    result[1][1] = s.input_space_size();
    for( BinWord a = 2; a < s.input_space_size(); a++)
    {
        result[a][0] = s.input_space_size();
        result[a][a] = s.input_space_size();
        for(BinWord x = 0, _max = s.input_space_size(), _msb = (BinWord)1 << cpp_msb(a) ; x < _max; x+=_msb)
        {
            std::vector<std::vector<BinWord> > xor_list(s.output_space_size(), std::vector<BinWord>(0));
            for(BinWord _ceil = x + _msb; x < _ceil; x++)
            {
                BinWord z = s[x] ^ s[x^a];
                for(auto x2 : xor_list[z]){
                    BinWord b = x ^ x2;
                    BinWord c = b ^ a;
                    if (c < a) {// Couldn't find a way to leverage this symetry :(
                        result[a][b] += 4;
                        result[a][c] += 4;
                        result[b][a] += 4;
                        result[b][c] += 4;
                        result[c][a] += 4;
                        result[c][b] += 4;
                    }
                }
                xor_list[z].push_back(x);
            }
        }
    }
    return result;
}
