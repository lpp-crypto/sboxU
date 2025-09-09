#include "boomerang.hpp"


// !SECTION! The BCT itself 


std::vector<Integer> cpp_bct_row(
    const cpp_S_box & s,
    const BinWord a)
{
    if (a == 0)
        return std::vector<Integer>(s.output_space_size(),s.input_space_size());

    std::vector<Integer> result(s.size(), 0);
    std::vector<std::vector<BinWord> > xor_list(s.output_space_size(), std::vector<BinWord>(0));
    result[0] = s.input_space_size();

    FOR_ENUMERATE_DIFFERENCE_COSETS(x,a,s.input_space_size())
    {
        BinWord z1 = s[x];
        BinWord z = z1 ^ s[x^a];

        for(auto z2 : xor_list[z]){
            result[ z1 ^ z2 ] += 4;
            if (z > 0)
                result[ z1 ^ z2 ^ z] += 4;
        }
        xor_list[z].push_back(z1);
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


// !SECTION! BCT Spectrum 


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
