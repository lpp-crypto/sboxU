/* Time-stamp: <2018-05-24 14:48:03 lperrin>
 *
 * LICENCE
 */

#include "sboxu_cpp_equiv.hpp"
using namespace boost::python;


// !SECTION! LEguess class


LEguess::LEguess(const unsigned int _target_size) :
    target_size(_target_size),
    latest_entries(),
    min_unset(1)
{
    partial_lut[0] = 0;
    is_set[0] = true;
}


LEguess::~LEguess()
{
    partial_lut.clear();
    is_set.clear();
}


std::vector<IOpair> LEguess::add_entry(const IOpair e)
{
    BinWord x = e.first, y = e.second;
    if ((x != 0) and (y == 0))
        throw ContradictionFound(x, y);
    else if ((x == 0) and (y != 0))
        throw ContradictionFound(x, y);
    else if (not is_set[x])
    {
        partial_lut[x] = y;
        is_set[x] = true;
        // propagating new value
        latest_entries.clear();
        std::map<BinWord, bool> previously_set(is_set);
        for (auto & entry : previously_set)
            if (entry.second)
            {
                BinWord
                    in_val = entry.first ^ x,
                    out_val = partial_lut[entry.first] ^ y;
                if (is_set[in_val])
                {

                    if (
                        (partial_lut[in_val] != out_val)
                        or
                        ((in_val != 0) and (out_val == 0))
                        )
                        throw ContradictionFound(in_val, out_val);
                }                    
                else
                {
                    partial_lut[in_val] = out_val;
                    is_set[in_val] = true;
                    latest_entries.push_back(IOpair(in_val, out_val));
                }
            }
        // updating the value of the minimal unset entry
        while ((min_unset < target_size) and is_set[min_unset])
            min_unset ++;
        return latest_entries;
    }
    else if (partial_lut[x] != y)
        throw ContradictionFound(x, y);
    return std::vector<IOpair>();
}


bool LEguess::is_entry_set(const BinWord x)
{
    return is_set[x];
}


unsigned int LEguess::min_u()
{
    return min_unset;
}


bool LEguess::complete()
{
    return (is_set.size() == target_size);
}


BinWord LEguess::img(const BinWord x)
{
    if (is_set[x])
        return partial_lut[x];
    else
        throw std::runtime_error("Trying to access entry which doesn't exist");
}


Sbox LEguess::lut()
{
    if (complete())
    {
        Sbox result(target_size, 0);
        for (auto & entry : partial_lut)
            result[entry.first] = entry.second;
        return result;
    }
    else
        throw std::runtime_error("Trying to return incomplete LUT");
}



// !SECTION! LEguessIterator class

LEguessIterator::LEguessIterator(const unsigned int _target_size,
                                 const std::vector<IOpair> _constraints) :
    target_size(_target_size),
    y(0),
    constraints(_constraints),
    base_guess(target_size),
    prepared_guess(target_size)
{
    unsigned long int two_pow_length = 1;
    guess_bit_length = 0;
    guess_mask = 0;
    while (two_pow_length < target_size)
    {
        guess_bit_length += 1;
        two_pow_length <<= 1;
        guess_mask = (guess_mask << 1) | 1;
    }
    if (guess_bit_length <= constraints.size())
        throw std::runtime_error("Guess too deep!");
    for (auto & c : constraints)
        base_guess.add_entry(c);
    x = base_guess.min_u();
}


bool LEguessIterator::prepare_successfully()
{
    y ++;
    if (y >= (target_size))
        return false;
    prepared_guess = LEguess(base_guess);
    try
    {
        prepared_guess.add_entry(IOpair(x, y));
    }
    catch (ContradictionFound & e)
    {
        return prepare_successfully();
    }
    return true;
}


LEguess LEguessIterator::get_prepared_guess()
{
    return prepared_guess;
}



LEguessIterator LEguessIterator::deeper_guess_gen()
{
    std::vector<IOpair> deeper_constraints(constraints);
    deeper_constraints.push_back(IOpair(x, y));
    return LEguessIterator(target_size, deeper_constraints);
}



// !SECTION! Actual algorithm LE

std::vector<Sbox> linear_equivalence_cpp(const Sbox f, const Sbox g)
{
    check_length(f);
    check_length(g);
    std::vector<Sbox> result;
    // if f(0) == 0 then it is necessary that g(0) == 0...
    if ((f[0] == 0) and (g[0] != 0))
        return result;
    // ... and vice-versa
    else if ((f[0] != 0) and (g[0] == 0))
        return result;
    // We look for A,B such that f = B o g o A
    Sbox f_inv = inverse(f), g_inv = inverse(g);
    std::vector<IOpair> constraints;
    if (f[0] != 0)
    {                           // in this case, we know that B(g[0]) = f[0]
        constraints.push_back(IOpair(g[0], f[0]));
    }
    std::vector<LEguessIterator> b_generators;
    b_generators.push_back(LEguessIterator(f.size(), constraints));
    while (true)
    {
        while ((b_generators.size() > 0) and (not b_generators.back().prepare_successfully()))
            b_generators.pop_back();
        if (b_generators.size() == 0)
                return std::vector<Sbox>(); // nothing to be found, we
                                            // have exhausted all
                                            // possible guesses
        LEguess
            a(f.size()),
            b = b_generators.back().get_prepared_guess();
        
        std::vector<IOpair> just_added, to_treat, next_to_treat;
        for (unsigned int x=0; x<f.size(); x++)
            if (b.is_entry_set(x))
                to_treat.push_back(IOpair(x, b.img(x)));
        bool
            conclusion_reached = false, // true if and only if we have
                                        // propagated the intial guess
                                        // as far as possible.
            contradiction_found = false; // true if and only if the
                                         // initial guess leads to
                                         // some inconsistency.
        // propagating the guess on B
        while (not conclusion_reached)
        {
            conclusion_reached = true;
            try
            {
                // B(y) = b_y => A(f_inv(b_y)) = g_inv(y)
                for (auto & entry : to_treat)
                {
                    BinWord y = entry.first, b_y = entry.second;
                    just_added = a.add_entry(IOpair(f_inv[b_y], g_inv[y]));
                    next_to_treat.insert(next_to_treat.end(),
                                         just_added.begin(),
                                         just_added.end());
                }
                to_treat = next_to_treat;
                next_to_treat.clear();
                if (to_treat.size() > 0)
                    conclusion_reached = false;
                // A(x) = a_x => B(g(a_x)) = f(x)
                for (auto & entry : to_treat)
                {
                    BinWord x = entry.first, a_x = entry.second;
                    just_added = b.add_entry(IOpair(g[a_x], f[x]));
                    next_to_treat.insert(next_to_treat.end(),
                                         just_added.begin(),
                                         just_added.end());
                }
                to_treat = next_to_treat;
                next_to_treat.clear();
                if (to_treat.size() > 0)
                    conclusion_reached = false;
            }
            catch (ContradictionFound & e)
            {
                conclusion_reached = true;
                contradiction_found = true;
            }
        }
        if (not contradiction_found)
        {
            if (not b.complete())
                // inconclusive: we will guess one more entry
                b_generators.push_back(b_generators.back().deeper_guess_gen());
            else
            {                   // we have everything we need
                std::vector<Sbox> result;
                Sbox A(f.size(), 0), B(b.lut());
                if (is_permutation_cpp(B))
                {
                    Sbox B_inv(inverse(B));
                    // A is rebuilt from scratch using that f = B o g o A
                    for (unsigned int x=0; x<B.size(); x++)
                        A[x] = g_inv[B_inv[f[x]]];
                    result.push_back(A);
                    result.push_back(B);
                    return result;
                }
                else
                    // inconclusive: we will guess one more entry
                    b_generators.push_back(b_generators.back().deeper_guess_gen());
                    
            }
        }
    }
}


// !SECTION! Representative LE class



LErepr::LErepr(const Sbox _f) :
    f(_f),
    f_inv(inverse(f)),
    contiguous_length(0),
    best_guess(f.size(), f.size()),
    a(0),
    b(0),
    a_inv(0),
    b_inv(0)
{
}


void LErepr::add_entry(const IOpair e)
{
    BinWord x = e.first, y = e.second;
    if (is_set[x])
    {
        if (partial_lut[x] != y)
            throw ContradictionFound(x, y);
    }
    if (is_inv_set[y])
    {
        if (partial_inv_lut[y] != x)
            throw ContradictionFound(x, y);
    }
    if ((not is_set[x]) or (not is_inv_set[y]))
    {
        partial_lut[x] = y;
        is_set[x] = true;
        partial_inv_lut[y] = x;
        is_inv_set[y] = true;
        if (is_set[contiguous_length + 1])
        {
            while ((contiguous_length < f.size()) and is_set[contiguous_length])
                contiguous_length ++;
            if (current_is_greater_than_best())
                throw ContradictionFound(0, 0);
        }
    }
}


Sbox LErepr::lut()
{
    return best_guess;
}


bool LErepr::current_is_greater_than_best()
{
    for (unsigned int x=0; x<contiguous_length; x++)
    {
        if (partial_lut[x] > best_guess[x])
            return true;
        else if (partial_lut[x] < best_guess[x])
            return false;
    }
    return false;
}


bool LErepr::update_state()
{    
    std::vector<IOpair> just_added;
    bool
        done_propagating = false,
        modified = false;
    while (not done_propagating)
    {
        done_propagating = true;
        // B(Rs(x)) = f(A(x))
        for (unsigned int x=0; x<f.size(); x++)
            if (a.is_entry_set(x) and is_set[x])
            {
                just_added = b.add_entry(IOpair(partial_lut[x],
                                                f[a.img(x)]));
                b_inv.add_entry(IOpair(f[a.img(x)],
                                       partial_lut[x]));
                if (just_added.size() > 0)
                    done_propagating = false;
            }
        // A(Rs^{-1}(y)) = f^{-1}(B(y))
        for (unsigned int y=0; y<f.size(); y++)
            if (b.is_entry_set(y) and is_inv_set[y])
            {
                just_added = a.add_entry(IOpair(partial_inv_lut[y],
                                                f_inv[b.img(y)]));
                a_inv.add_entry(IOpair(f_inv[b.img(y)],
                                       partial_inv_lut[y]));
                if (just_added.size() > 0)
                    done_propagating = false;
            }
        // Rs(x) = B^{-1}(f(A(x)))
        for (unsigned int x=0; x<f.size(); x++)
            if (a.is_entry_set(x)
                and
                b_inv.is_entry_set(f[a.img(x)]))
            {
                unsigned int y = b_inv.img(f[a.img(x)]);
                if (not is_set[x])
                    done_propagating = false;
                add_entry(IOpair(x, y));
            }
        if (not done_propagating)
            modified = true;
    }
    return modified;
}


unsigned int LErepr::min_n_a()
{
    // D_A is the set of points x such that A(x) is known but not
    // Rs(x)
    for (auto & entry : a.is_set)
        if (entry.second and not is_set[entry.first])
            return entry.first;
    return f.size();
}


unsigned int LErepr::min_n_b()
{
    // D_A is the set of points x such that B(x) is known but
    // not Rs^{-1}(x)
    for (auto & entry : b.is_set)
        if (entry.second and not is_inv_set[entry.first])
            return entry.first;
    return f.size();
}


unsigned int LErepr::min_u_a()
{
    // U_A is the set of points x such that A(x) is not known
    return a.min_u();
}


unsigned int LErepr::min_u_b()
{
    // U_A is the set of points x such that B(x) is not known
    return b.min_u();
}


void LErepr::initialize()
{
    // IMPORTANT REMARK:
    // 
    // We look for Rs which is linearly equivalent to f and the linear
    // permutations A and B such that Rs = B^{-1} o f o A. These
    // notations are those used in the paper of Biryukov et al but
    // they are different from those in are_linear_equivalent_cpp
    // (although they are also borrowed from the same paper).

    std::vector<IOpair> constraints_a;
    if (f[0] != 0)
        constraints_a.push_back(IOpair(1, f_inv[0]));
    std::vector<LEguessIterator> a_generators;
    a_generators.push_back(LEguessIterator(f.size(),
                                           constraints_a));
    // In this loop, all guesses for A are explored. The guesses are
    // built iteratively: we start by guessing only one entry. If it
    // leads to a contradiction, we move on to the next value. If no
    // conclusion is reached, then we add an additional variable to
    // guess. We effectively explore a tree corresponding to all
    // possible A by making a decision on each branch as early as
    // possible.
    while (true)
    {
        while ((a_generators.size() > 0) and (not a_generators.back().prepare_successfully()))
            a_generators.pop_back();
        if (a_generators.size() == 0)
            return ; // nothing more to be found, we have exhausted
                     // all possible guesses.
        // reinitializing guesses
        a = a_generators.back().get_prepared_guess();
        b = LEguess(f.size());
        a_inv = LEguess(f.size());
        b_inv = LEguess(f.size());
        // reinitialize state
        partial_lut.clear();
        partial_inv_lut.clear();
        is_set.clear();
        is_inv_set.clear();
        // adding first two entries to Rs and updating a_inv
        bool contradiction_found = false;
        try
        {
            if (f[0] == 0)
            {
                add_entry(IOpair(0, 0));
                add_entry(IOpair(1, 1));
            }
            else
            {
                add_entry(IOpair(0, 1));
                add_entry(IOpair(1, 0));
            }
            for (auto & entry : a.is_set)
                if (entry.second)
                    a_inv.add_entry(IOpair(a.img(entry.first), entry.first));
            (void)update_state();
        }
        catch (ContradictionFound & e) {contradiction_found = true;}
        // Main loop
        unsigned int x = 0, y = 0;
        while ((min_n_a() < f.size())     // N_A is not empty
               and (not contradiction_found))
        {
            x = min_n_a();
            y = min_u_b();
            try
            {
                add_entry(IOpair(x, y));
                b.add_entry(IOpair(y, f[a.img(x)]));
                b_inv.add_entry(IOpair(f[a.img(x)], y));
                (void)update_state();
            }
            catch (ContradictionFound & e) {contradiction_found = true; }
            while ((min_n_a() == f.size()) // N_A is empty
                   and (min_n_b() < f.size()) // N_B is not empty
                   and (not contradiction_found))
            {
                try
                {
                    x = min_u_a();
                    y = min_n_b();
                    add_entry(IOpair(x, y));
                    a.add_entry(IOpair(x, f_inv[b.img(y)]));
                    a_inv.add_entry(IOpair(f_inv[b.img(y)], x));
                    (void)update_state();
                }
                catch (ContradictionFound & e) {contradiction_found = true; }
            }
        }
        if (not contradiction_found) // then we have found a new candidate for Rs!
        {

            bool is_new_candidate_valid = true;
            // testing if it is fully specified
            for (unsigned int x=0; x<f.size(); x++)
                if (not is_set[x])
                {
                    // if Rs is not fully specified, we need to guess deeper
                    a_generators.push_back(a_generators.back().deeper_guess_gen());
                    is_new_candidate_valid = false;
                    break;
                }
            // if it is fully specified and smaller, replacing best candidate
            if (is_new_candidate_valid and (not current_is_greater_than_best()))
            {
                for (unsigned int x=0; x<f.size(); x++)
                    best_guess[x] = partial_lut[x];
            }
        }
    }
}

Sbox le_class_representative_cpp(const Sbox f)
{
    check_length(f);
    list result;
    LErepr repr(f);
    repr.initialize();
    return repr.lut();
}



// !SECTION! Vector extraction


std::vector<BinWord> extract_vector_cpp(
    const std::vector<BinWord> z,
    const BinWord a)
{
    // building indicator function
    std::map<BinWord, bool> indicator;
    for(auto & x : z)
    {
        indicator[x] = true;
    }
    std::vector<BinWord> result;
    for(auto & x : z)
    {
        long int y = x ^ a;
        if ((x < y) and (indicator.find(y) != indicator.end()))
            result.push_back(x);
    }
    return result;
}


std::vector<BinWord> super_extract_vector_cpp(
    const std::vector<BinWord> z,
    const BinWord a)
{
    std::set<BinWord> indicator(z.begin(), z.end());
    std::vector<BinWord> result; //(1, 0);
    result.reserve(1 + (z.size()/2)) ;
    for(auto & x : z)
        if (a < x)
        {
            BinWord y = x ^ a;
            if (x < y)
            {
                if (indicator.find(y) != indicator.end())
                    result.push_back(x);
            }
        }
    indicator.clear();
    return result;
}


std::vector<std::vector<BinWord> > extract_bases_rec(
    const std::vector<BinWord> base,
    const std::vector<BinWord> z,
    const Integer dimension,
    const Integer word_length)
{
    std::vector<std::vector<BinWord> > result, tmp;
    std::vector<BinWord> new_base, z_a;
    bool branch_is_exhausted = false;

    if (z.size() == 0)
    {   // branch is done if z is empty
        branch_is_exhausted = true;
    }
    else if (base.size() < dimension)
    {   // branch is done if z does not contain enough element to form
        // a base of dimension `dimension`
        branch_is_exhausted = (z.size() < ((1 << (dimension - base.size())) - 1)) ;
    }
    if (branch_is_exhausted && (base.size() >= dimension))
    {   // if the search stopped but the base is large enough, we keep it
        result.push_back(base);
    }
    else
    {   // otherwise, the search continues
        BinWord mask = 0;
        for (unsigned int i=word_length-dimension+base.size()+1; i<word_length; i++)
            mask |= (1 << i) ;
        for (unsigned int i=0; i<z.size(); i++)  // recursively extracting
        { 
            BinWord a = z[i];
            if ((a > base.back()) && ((a & mask) == 0))
            {
                new_base.assign(base.begin(), base.end());
                new_base.push_back(a);
                z_a = super_extract_vector_cpp(std::vector<BinWord>(z.begin()+i+1, z.end()), a);
                tmp = extract_bases_rec(new_base, z_a, dimension, word_length);
                result.insert(result.end(), tmp.begin(), tmp.end());
                tmp.clear();
                z_a.clear();
            }
        }
    }
    return result;
}


void extract_bases_starting_with(
    const std::vector<BinWord> starting_vectors,
    std::vector<std::vector<BinWord> > & result,
    const std::vector<BinWord> z,
    const Integer dimension,
    const Integer word_length)
{    
    std::vector<std::vector<BinWord> > tmp, valid_starts;
    std::vector<BinWord> new_base, z_a;

    // looping over starting points
    for (unsigned int i=0; i<starting_vectors.size(); i++)
    {
        BinWord a = starting_vectors[i] ;
        // if ((a & mask) == 0)
        {
            z_a = super_extract_vector_cpp(std::vector<BinWord>(z.begin()+i+1, z.end()), a);
            new_base.assign(1, a);
            tmp = extract_bases_rec(new_base, z_a, dimension, word_length);
            result.insert(result.end(), tmp.begin(), tmp.end());
            tmp.clear();
            z_a.clear();
        }
    }
}


std::vector<std::vector<BinWord> > extract_bases_cpp(
    std::vector<BinWord> z,
    const Integer dimension,
    const Integer word_length,
    const Integer n_threads)
{
    std::vector<std::vector<BinWord> > result;
    std::vector<std::thread> threads;
    std::vector<std::vector<std::vector<BinWord> > > local_results(
        n_threads,
        std::vector<std::vector<BinWord> >());
    std::vector<std::vector<BinWord> > all_starting_vectors(n_threads, std::vector<BinWord>());

    // ensuring that the list is initially sorted
    std::sort(z.begin(), z.end());
    
    // placing all relevant vectors in different buckets
    BinWord mask = 0;
    for (unsigned int i=word_length-dimension+1; i<word_length; i++)
        mask |= (1 << i);
    unsigned int counter = 0; 
    for (auto &a : z)
        if ((a != 0) and ((a & mask) == 0))
        {
            all_starting_vectors[counter % n_threads].push_back(a);
            counter ++;
        }
    // assigning each bucket to a different thread
    for (unsigned int i=0; i<n_threads; i++)
    {
        threads.push_back(std::thread(extract_bases_starting_with,
                                      all_starting_vectors[i],
                                      std::ref(local_results[i]),
                                      z,
                                      dimension,
                                      word_length));
    }
    // regrouping results
    for (unsigned int i=0; i<n_threads; i++)
    {
        threads[i].join();
        result.insert(result.end(),
                      local_results[i].begin(),
                      local_results[i].end());
    }    
    return result;
}

