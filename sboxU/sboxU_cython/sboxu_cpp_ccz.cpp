/* Time-stamp: <2021-08-12 17:16:07 lperrin>
 *
 * LICENCE
 */

#include "sboxu_cpp_ccz.hpp"


#define VECTOR_EXTRACT_ALL_DIMS 0
#define VECTOR_EXTRACT_FIXED_DIM 1
#define VECTOR_EXTRACT_JUST_ONE 2

#define VECTOR_EXTRACT_DEFAULT VECTOR_EXTRACT_ALL_DIMS


// !SECTION! Linear extraction


std::vector<BinWord> extract_vector_cpp(
    const std::vector<BinWord> & z,
    const BinWord a)
{
    // building indicator function
    std::map<BinWord, bool> indicator;
    for(auto & x : z)
        if (x > a)
        {
            indicator[x] = true;
        }
    std::vector<BinWord> result;
    for(auto & x : z)
        if (x > a)
        {
            BinWord y = x ^ a;
            if ((x < y) and (indicator.find(y) != indicator.end()))
                result.push_back(x);
        }
    return result;
}




class MSBSpectrum
{
private:
    std::vector<BinWord> clz_switching_positions ;
    Integer dimension ;
        
public:
    MSBSpectrum(const std::vector<BinWord> & z, const Integer _dimension) :
        dimension(_dimension)
    {
        clz_switching_positions.push_back(0);
        unsigned int previous_switch_clz = __builtin_clz(z[0]);
        for (unsigned int i=1; i<z.size(); i++)
        {
            if (previous_switch_clz > __builtin_clz(z[i]))
            {
                // if the number of leading zeroes in the last entry
                // in clz_switching_positions is lower than in z[i],
                // then i is a position at which the number of leading
                // zeroes switches
                clz_switching_positions.push_back(i);
                previous_switch_clz = __builtin_clz(z[i]);
            }
        }
        clz_switching_positions.push_back(z.size());
    };

    unsigned int class_size(const unsigned int i) const
    {
        return clz_switching_positions[i+1] - clz_switching_positions[i] ;
    }    

    unsigned int slice_start(const unsigned int i) const
    {
        return clz_switching_positions[i] ;
    }
    
    unsigned int slice_end(const unsigned int i) const
    {
        return clz_switching_positions[i+1] ;
    }

    unsigned int size() const
    {
        return clz_switching_positions.size() - 1;
    }

    bool contains_MSB_sequence_starting_from(const unsigned int i) const
    {
        if (dimension <= 0)
            return true;
        else
        {
            unsigned int j = 0;
            for (unsigned int k=i; k<size(); k++)
            {
                if (class_size(k) >= (((unsigned int)1) << j))
                    j += 1;
                if ((dimension - j) > (size() - k))
                    return false;
                if (j >= (unsigned int)dimension)
                    return true;
            }
            return false;
        }
    }
};



std::vector<BinWord> super_extract_vector_cpp(
    const std::vector<BinWord> &z,
    const MSBSpectrum & spec,
    const BinWord a,
    const unsigned int slice_of_a)
{
    std::vector<BinWord> result;
    result.reserve((z.size() - spec.slice_end(slice_of_a))/2) ;
    for (unsigned int slice_index=slice_of_a+1; slice_index<spec.size(); slice_index++)
    {
        for (unsigned int i=spec.slice_start(slice_index); i<spec.slice_end(slice_index); i++)
        {
            const BinWord x = z[i];
            const BinWord y = x ^ a;
            if (x < y)
            {
                if (std::binary_search(z.begin()+i+1,
                                       z.begin()+spec.slice_end(slice_index),
                                       y))
                    result.push_back(x) ;
            }            
        }
    }
    return result;
}



std::vector<std::vector<BinWord> > extract_bases_rec(
    const std::vector<BinWord> &z,
    const Integer dimension,
    const Integer end_code)
{
    std::vector<std::vector<BinWord> > result;

    if (z.size() == 0)
    {
        // if the word list is empty then there cannot be any vector
        // space in it. However, if we were so deep in the recursion
        // that we went beyond the wished dimension, then we add an
        // empty vector to `result` so that this vector can receive
        // the basis vector in the previous recursive calls.
        if ((end_code == VECTOR_EXTRACT_ALL_DIMS) && (dimension <= 0))
        {
            result.emplace_back() ;
        }
        return result;
    }
    else if ((dimension == 1) && (end_code != VECTOR_EXTRACT_ALL_DIMS))
    {
        // if the list is not empty and if we are almost done (i.e. we
        // only want spaces of dimension 1)
        if (end_code == VECTOR_EXTRACT_JUST_ONE)
        {
            // since we only care about finding 1 basis, we only
            // return the first one.
            result.emplace_back(1, z[0]);
        }
        else if (end_code == VECTOR_EXTRACT_FIXED_DIM)
        {
            // since we want all the basis with given dimension, we
            // have to return all the elements in z.
            result.reserve(z.size());
            for (auto & a: z)
                result.emplace_back(1, a);
        }
        return result;
    }
    else
    {
        // at that point, z is not empty and either the dimension is
        // >1 or we are extracting all dimensions.
        
        // If the set does not satisfy the MSB criteria, we stop
        const MSBSpectrum spec(z, dimension) ;
        if (spec.size() < (unsigned int)dimension)
        {
            if ((end_code == VECTOR_EXTRACT_ALL_DIMS) && (dimension <= 0))
            {
                // if we were going deeper in the recursion to find
                // spaces of dimension greater than the one sought
                // then we let all current list elements be the end of
                // GJBs.
                result.reserve(z.size());
                for (auto & a: z)
                    result.emplace_back(1, a);
            }
            return result;
        }
         

        // Otherwise, the search continues. We loop recursively over
        // all slices independently, a slice corresponding to the set
        // of all elements having the same MSB
        const unsigned int min_size_extracted = (1 << (dimension-1)) - 1;
        for (unsigned int slice_index=0; slice_index<spec.size(); slice_index++)
        {
            unsigned int
                start = spec.slice_start(slice_index),
                end   = spec.slice_end(slice_index);
            if (spec.contains_MSB_sequence_starting_from(slice_index))
            {                
                for (unsigned int i=start; i<end; i++)
                {
                    const BinWord a = z[i];
                    std::vector<BinWord> z_a = std::move(super_extract_vector_cpp(
                                                                   z,
                                                                   spec,
                                                                   a,
                                                                   slice_index));
                    if (z_a.size() >= min_size_extracted)
                    {
                        std::vector<std::vector<BinWord> > tmp = std::move(
                            extract_bases_rec(z_a,
                                              dimension-1,
                                              end_code)
                            );
                        result.reserve(result.size() + tmp.size()) ;
                        for (auto & base : tmp)
                        {
                            std::vector<BinWord> new_base(1, a);
                            new_base.insert(new_base.end(), base.begin(), base.end());
                            result.push_back(new_base);
                            if (end_code == VECTOR_EXTRACT_JUST_ONE)
                                return result;
                        }
                    }
                    else if ((end_code == VECTOR_EXTRACT_ALL_DIMS) && (dimension <= 0))
                    {
                        // if we were going deeper in the recursion to
                        // find spaces of dimension greater than the
                        // one sought then we let the current list
                        // element be the end of a GJB.
                        result.emplace_back(1, a);
                    }
                }
            }
            else if ((end_code == VECTOR_EXTRACT_ALL_DIMS) && (dimension <= 0))
            {
                // if we were going deeper in the recursion to find
                // spaces of dimension greater than the one sought
                // then we let all current slide elements be the end of
                // GJBs.
                for (unsigned int i=start; i<end; i++)
                    result.emplace_back(1, z[i]);
            }
        }
    }
    return result;
}


void extract_bases_starting_with(
    const std::vector<BinWord> & starting_vectors,
    std::vector<std::vector<BinWord> > & result,
    const std::vector<BinWord> & z,
    const MSBSpectrum & spec,
    const Integer dimension,
    const Integer end_code)
{    
    BinWord min_size_extracted = (1 << (dimension - 1)) - 1 ;
    unsigned int
        index_in_z = 0,
        slice_index = 0;
    // precomputing which slices are useful
    std::vector<bool> slice_is_good(spec.size(), false);
    for (unsigned int slice_index=0; slice_index<spec.size(); slice_index++)
        slice_is_good[slice_index] = spec.contains_MSB_sequence_starting_from(slice_index) ;

    
    // looping over starting points
    for (unsigned int i=0; i<starting_vectors.size(); i++)
    {
        BinWord a = starting_vectors[i] ;
        while (z[index_in_z] != a)
            index_in_z ++ ;
        while (index_in_z >= spec.slice_end(slice_index))
            slice_index ++;
        if (slice_is_good[slice_index])
        {        
            // Vector Extraction
            std::vector<BinWord> z_a = std::move(super_extract_vector_cpp(
                                                     z,
                                                     spec,
                                                     a,
                                                     slice_index));
            // Continuing if the result of the extraction is big enough
            if (z_a.size() >= min_size_extracted)
            {
                std::vector<std::vector<BinWord> > tmp =
                    std::move(extract_bases_rec(z_a, dimension-1, end_code));
                result.reserve(result.size() + tmp.size()) ;
                for (auto & base : tmp)
                {
                    std::vector<BinWord> new_base(1, a);
                    new_base.insert(new_base.end(), base.begin(), base.end());
                    result.push_back(new_base);
                    if (end_code == VECTOR_EXTRACT_JUST_ONE)
                        return ;                
                }
            }
        }
        else
        {
            break;
        }
    }
}


std::vector<std::vector<BinWord> > extract_bases_cpp(
    std::vector<BinWord> & z,
    const Integer dimension,
    const Integer word_length,
    Integer n_threads,
    const std::string end_condition)
{
    std::vector<std::vector<BinWord> > result;
    if (z.size() == 0)
        return result;
    
    Integer end_code = VECTOR_EXTRACT_DEFAULT;
    if (end_condition.compare("all dimensions") == 0)
        end_code = VECTOR_EXTRACT_ALL_DIMS;
    else if (end_condition.compare("fixed dimension") == 0)
        end_code = VECTOR_EXTRACT_FIXED_DIM;
    else if (end_condition.compare("just one") == 0)
    {
        n_threads = 1;
        end_code = VECTOR_EXTRACT_JUST_ONE;
    }
    
    std::vector<std::thread> threads;
    std::vector<std::vector<std::vector<BinWord> > > local_results(
        n_threads,
        std::vector<std::vector<BinWord> >());
    std::vector<std::vector<BinWord> > all_starting_vectors(n_threads, std::vector<BinWord>());

    // ensuring that the list is initially sorted
    std::sort(z.begin(), z.end());

    // MSB spectrum of z
    const MSBSpectrum spec(z, dimension);
    
    // placing all relevant vectors in different buckets
    for (unsigned int i=0; i<all_starting_vectors.size(); i++)
        all_starting_vectors[i].reserve(z.size() / n_threads);
    unsigned int bucket_index = 0;
    for (unsigned int slice_index=0; slice_index<spec.size(); slice_index++)
    {
        if (spec.contains_MSB_sequence_starting_from(slice_index))
        {
            for (unsigned int i=spec.slice_start(slice_index); i<spec.slice_end(slice_index); i++)
            {
                BinWord a = z[i] ;
                if (a != 0)
                {
                    all_starting_vectors[bucket_index].push_back(a);
                    bucket_index = (bucket_index+1) % n_threads;
                }
            }
        }
        else
        {
            break;
        }
    }

    // assigning each bucket to a different thread
    for (int i=0; i<n_threads; i++)
    {
        threads.emplace_back(extract_bases_starting_with,
                             all_starting_vectors[i],
                             std::ref(local_results[i]),
                             z,
                             spec,
                             dimension,
                             end_code);
    }
    // regrouping results
    for (int i=0; i<n_threads; i++)
    {
        threads[i].join();
        result.insert(result.end(),
                      local_results[i].begin(),
                      local_results[i].end());
    }    
    return result;
}


// !SECTION! Affine Extraction

std::vector<BinWord> affine_pre_process(
    const BinWord c,
    const std::vector<BinWord> & z)
{
    std::vector<BinWord> result;
    result.reserve(z.size());
    for (auto & x : z)
    {
        if (x > c)
            result.push_back(x ^ c) ;
    }
    // ensuring that the list is initially sorted
    std::sort(result.begin(), result.end());
    return result;
}


void extract_affine_bases_starting_with(
    const std::vector<BinWord> & starting_vectors,
    std::vector<std::vector<BinWord> > & result,
    const std::vector<BinWord> & z,
    const Integer dimension,
    const Integer end_code)
{    
    // looping over starting points
    for (auto &a : starting_vectors)
    {
        // pre-processing
        const std::vector<BinWord> pre_processed_z = affine_pre_process(a, z) ;
        // linear extraction
        std::vector<std::vector<BinWord> > tmp =
            std::move(extract_bases_rec(pre_processed_z, dimension, end_code));
        result.reserve(result.size() + tmp.size()) ;
        for (auto & base : tmp)
        {
            std::vector<BinWord> new_base(1, a);
            new_base.insert(new_base.end(), base.begin(), base.end());
            result.push_back(new_base);
            if (end_code == VECTOR_EXTRACT_JUST_ONE)
                return ;
        }
    }
}



std::vector<std::vector<BinWord> > extract_affine_bases_cpp(
    std::vector<BinWord> & z,
    const Integer dimension,
    const Integer word_length,
    Integer n_threads,
    const std::string end_condition)
{
    std::vector<std::vector<BinWord> > result;
    if (z.size() == 0)
        return result;
        
    Integer end_code = VECTOR_EXTRACT_DEFAULT;
    if (end_condition.compare("all dimensions") == 0)
        end_code = VECTOR_EXTRACT_ALL_DIMS;
    else if (end_condition.compare("fixed dimension") == 0)
        end_code = VECTOR_EXTRACT_FIXED_DIM;
    else if (end_condition.compare("just one") == 0)
    {
        n_threads = 1;
        end_code = VECTOR_EXTRACT_JUST_ONE;
    }
    
    std::vector<std::thread> threads;
    std::vector<std::vector<std::vector<BinWord> > > local_results(
        n_threads,
        std::vector<std::vector<BinWord> >());
    std::vector<std::vector<BinWord> > all_starting_vectors(n_threads, std::vector<BinWord>());

    std::sort(z.begin(), z.end());

    // placing all relevant vectors in different buckets
    for (unsigned int i=0; i<all_starting_vectors.size(); i++)
        all_starting_vectors[i].reserve(z.size() / n_threads);
    unsigned int bucket_index = 0;
    for (auto &a : z)
    {
        all_starting_vectors[bucket_index].push_back(a);
        bucket_index = (bucket_index+1) % n_threads;
    }

    // assigning each bucket to a different thread
    for (int i=0; i<n_threads; i++)
    {
        threads.emplace_back(extract_affine_bases_starting_with,
                             all_starting_vectors[i],
                             std::ref(local_results[i]),
                             z,
                             dimension,
                             end_code);
    }
    // regrouping results
    for (int i=0; i<n_threads; i++)
    {
        threads[i].join();
        result.insert(result.end(),
                      local_results[i].begin(),
                      local_results[i].end());
    }    
    return result;
}



// !SECTION! Sigma multiplicities



void sigma_multiplicities_local_cpp(
    std::map<BinWord, Integer> & result,
    const Sbox f,
    const Integer k,
    const Integer n,
    const unsigned int s_min,
    const unsigned int s_max)
{
    for (unsigned int s=s_min; s<s_max; s++)
    {
        Integer multiplicity = 0;
        for (unsigned int b=0; b<f.size(); b++)
        {
            int sign = (scal_prod_cpp(b, s) == 0) ? 1 : -1;
            std::vector<Integer> w = walsh_spectrum_coord(component_cpp(b, f));
            for (auto &l : w)
            {
                Integer pow = l;
                for (unsigned int i=1; i<k; i++)
                    pow *= l;
                multiplicity += sign * pow; 
            }
        }
        result[multiplicity >> (2*n)] += 1;
    }
}


std::map<BinWord, Integer> sigma_multiplicities_cpp(const Sbox f,
                                                    const Integer k,
                                                    const Integer n,
                                                    const Integer n_threads)
{
    check_length_cpp(f);
    std::map<BinWord, Integer> result;
    if (n_threads == 1)
    {
        // small S-Box
        sigma_multiplicities_local_cpp(std::ref(result), f, k, n, 0, f.size());
    }
    else
    {
        std::vector<std::thread> threads;
        std::vector<std::map<BinWord, Integer> > local_counts(n_threads);
        unsigned int slice_size = f.size()/n_threads;
        for (unsigned int i=0; i<n_threads; i++)
        {
            unsigned int
                lower_bound = i*slice_size,
                upper_bound = (i+1)*slice_size;
            if (upper_bound > f.size())
                upper_bound = f.size();
            threads.push_back(std::thread(sigma_multiplicities_local_cpp,
                                          std::ref(local_counts[i]),
                                          f,
                                          k,
                                          n,
                                          lower_bound,
                                          upper_bound));

        }
        for (unsigned int i=0; i<n_threads; i++)
        {
            threads[i].join();
            for (auto &entry : local_counts[i])
                result[entry.first] += entry.second ;
        }
    }
    return result;
}


