#include "invariants.hpp"

// !SECTION! Ortho-derivative


// !SUBSECTION! The ortho-derivative itself

cpp_S_box cpp_ortho_derivative(const cpp_S_box s)
{
    std::vector<BinWord> result(s.input_space_size(), 0);
    for (unsigned int a=1; a<s.input_space_size(); a++)
    {
        // getting the hyperplane
        std::vector<Integer> row(cpp_ddt_row(s, a));
        std::vector<Integer> hyperplane;
        hyperplane.reserve(s.size() / 2);
        for(unsigned int b=1; b<s.input_space_size(); b++)
            // !TODO! use the hyperplane iterator instead? 
            if (row[b] != row[0])
                hyperplane.push_back(b);

        // we return an empty list if the function is not APN, which
        // is equivalent to all rows having exactly half of their
        // elements be non-zero
        if (hyperplane.size() < (s.input_space_size()/2))
            return cpp_empty_S_box();

        // bruteforcing "ortho" until it is orthogonal to all elements
        // in the hyperplane
        BinWord ortho = 1;
        bool found = false;
        while ((not found) and (ortho < s.input_space_size()))
        {
            found = true;
            for(auto &b : hyperplane)
                if (cpp_scal_prod(ortho, b) == 0)
                {
                    found = false;
                    break;
                }
            if (not found)
                ortho += 1;
        }
        // if we couldn't find an element orthogonal to the
        // hyperplane, then it is not a hyperplane and we return an
        // empty function
        if (found)
            result[a] = ortho;
        else
            return cpp_empty_S_box();
    }
    return cpp_S_box(result);
}



// !SECTION! Sigma multiplicities



void cpp_sigma_multiplicities_local(
    cpp_Spectrum & result,
    const cpp_S_box f,
    const Integer k,
    const unsigned int s_min,
    const unsigned int s_max)
{
    for (unsigned int s=s_min; s<s_max; s++)
    {
        Integer multiplicity = 0;
        for (unsigned int b=0; b<f.output_space_size(); b++)
        {
            int sign = (cpp_scal_prod(b, s) == 0) ? 1 : -1;
            std::vector<Integer> w = cpp_walsh_transform(f.component(b));
            for (auto &l : w)
            {
                Integer pow = l;
                for (unsigned int i=1; i<k; i++)
                    pow *= l;
                multiplicity += sign * pow; 
            }
        }
        result.incr(multiplicity >> (2*f.get_input_length()));
    }
}


cpp_Spectrum cpp_sigma_multiplicities(
    const cpp_S_box f,
    const Integer k,
    const Integer n_threads)
{
    cpp_Spectrum result;
    if (n_threads == 1)         // single thread
    {
        cpp_sigma_multiplicities_local(std::ref(result), f, k, 0, f.size());
    }
    else                        // strict multi-threading
    {
        std::vector<std::thread> threads;
        std::vector<cpp_Spectrum> local_counts(n_threads);
        unsigned int lower_bound = 0;
        for (unsigned int i=0; i<n_threads; i++)
        {
            // Will break on 32-bit arch if nthreads*s.size >= 1 << 32
            BinWord upper_bound = ((i+1)*f.output_space_size())/n_threads;
            threads.push_back(std::thread(cpp_sigma_multiplicities_local,
                                          std::ref(local_counts[i]),
                                          f,
                                          k,
                                          lower_bound,
                                          upper_bound));

        }
        for (unsigned int i=0; i<n_threads; i++)
        {
            threads[i].join();
            result += local_counts[i];
        }
    }
    return result;
}
