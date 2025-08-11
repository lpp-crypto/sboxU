#include "boomerang.hpp"



// !SECTION! The BCT itself 


std::vector<Integer> bct_row(
    const cpp_S_box & s,
    const cpp_S_box & s_inv,
    const BinWord a)
{
    if (a == 0)
        return std::vector<Integer>(s.output_space_size(),s.input_space_size());

    std::vector<Integer> result(s.size(), 0);
    std::vector<std::vector<BinWord> > xor_list(s.output_space_size(), std::vector<BinWord>(0));

    result[0] = s.input_space_size();
    BinWord max_bit = s.input_space_size()/2 ;
    for (unsigned int i=0; i<max_bit; i++)
    {
        BinWord x = i;
        BinWord y = x^a;
        if (y < x){
            x+=max_bit;
            y+=max_bit;
        }
        BinWord z = s[x] ^ s[y];
        xor_list[z].push_back(s[x]);
        result[z] += 2;
    }


    for (BinWord z=0; z<s.size(); z++)
    {
        auto collision_list = xor_list[z];
        size_t len = collision_list.size();
        for(size_t i = 0; i < len; i++)
        {
            BinWord z1 = collision_list[i];
            for(size_t j = i+1; j < len; j++)
            {
                BinWord z2 = collision_list[j];
                result[ z1 ^ z2 ] += 4;
                if( z > 0 )
                    result[ z1 ^ z2 ^ z ] += 4;
            }
        }
    }
    return result;
}


std::vector< std::vector<Integer>> cpp_bct(const cpp_S_box & s)
{
	std::vector< std::vector<Integer>> table_bct ;
        table_bct.reserve(s.input_space_size());
	cpp_S_box s_inv = s.inverse() ;
	for(unsigned int a = 0; a < s.input_space_size(); a++)
		table_bct.push_back(bct_row(s,s_inv,a)) ;
	return table_bct ;
}


// !SECTION! BCT Spectrum 

void bct_rows_count(
    cpp_Spectrum &result,
    const cpp_S_box & s,
    const cpp_S_box & s_inv,
    const BinWord a_min,
    const BinWord a_max)
{
    for (unsigned int a=a_min; a<a_max; a++)
    {
        std::vector<Integer> row = bct_row(s, s_inv, a);
        for (unsigned int i=1; i<row.size(); i++) // we start at 1, not 0
            result.incr(row[i]);
    }
}


cpp_Spectrum cpp_boomerang_spectrum(const cpp_S_box & s, const unsigned int n_threads)
{
    cpp_S_box s_inv = s.inverse() ;
    cpp_Spectrum count;
    if (n_threads == 1)         // single thread
    {
        bct_rows_count(std::ref(count), s, s_inv, 1, s.size());
    }
    else                        // multiple threads
    {
        std::vector<std::thread> threads;
        std::vector<cpp_Spectrum> local_counts(n_threads);
        BinWord lower_bound = 1;
        for (unsigned int i=0; i<n_threads; i++)
        {
            // Will break on 32-bit arch if nthreads*s.size >= 1 << 32
            BinWord upper_bound = ((i+1)*s.output_space_size())/n_threads;
            threads.push_back(std::thread(bct_rows_count,
                                          std::ref(local_counts[i]),
                                          s,
                                          s_inv,
                                          lower_bound,
                                          upper_bound));

            lower_bound = upper_bound;
        }
        for (unsigned int i=0; i<n_threads; i++)
        {
            threads[i].join();
            count += local_counts[i];
        }
    }
    return count;
}
