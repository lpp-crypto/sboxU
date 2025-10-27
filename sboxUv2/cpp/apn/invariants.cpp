#include "invariants.hpp"

// !SECTION! Invariant computation


// !SUBSECTION! The ortho-derivative itself

cpp_S_box cpp_ortho_derivative(const cpp_S_box &s)
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



// !SUBSECTION! Sigma multiplicities



void cpp_sigma_multiplicities_local(
    cpp_Spectrum & result,
    const cpp_S_box &f,
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
    const cpp_S_box &f,
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
            lower_bound = upper_bound;

        }
        for (unsigned int i=0; i<n_threads; i++)
        {
            threads[i].join();
            result += local_counts[i];
        }
    }
    return result;
}



// !SECTION! Putting invariants together


std::string cpp_apn_ea_mugshot(
    const cpp_Spectrum &abs_walsh_spec,
    const cpp_Spectrum &deg_spec,
    const cpp_Spectrum &sig_mult,
    const cpp_Spectrum &thk_spec
    )
{
    if (deg_spec.maximum() == 2)
        throw std::runtime_error("cpp_apn_ea_mughot(...spectra...) shouldn't be called on a quadratic function");
    std::stringstream result;
    result << abs_walsh_spec.content_string_repr()
           <<       deg_spec.content_string_repr()
           <<       sig_mult.content_string_repr()
           <<       thk_spec.content_string_repr();
    return result.str();
}


std::string cpp_apn_ea_mugshot(
    const cpp_S_box &s,
    const unsigned int n_threads
    )
{
    if (cpp_algebraic_degree(s) == 2)
    {
        // case of a quadratic functions: we rely on the ortho-derivative
        cpp_S_box o = cpp_ortho_derivative(s);
        // the linear equivalence class of the ortho-derivative is an
        // EA-equivalence class invariant for quadratic APN
        // functions. Thus, we use the (signed) Walsh spectrum and the
        // differential spectrum of the ortho-derivative. We don't use
        // its degree spectrum as it is always the same.
        std::stringstream result;
        result <<        cpp_walsh_spectrum(o, n_threads).content_string_repr()
               << cpp_differential_spectrum(o, n_threads).content_string_repr();
        return result.str();
    }
    else
        return cpp_apn_ea_mugshot(
            cpp_absolute_walsh_spectrum(s, n_threads),
            cpp_degree_spectrum(s),
            cpp_sigma_multiplicities(s, 4, n_threads),
            cpp_thickness_spectrum(s, n_threads)
            );
}
