/* Time-stamp: <2023-01-04 16:11:57 lperrin>
 *
 * LICENCE
 */

#include "sboxu_cpp_equiv_approx.hpp"


// !SECTION! LEguess class


LEguessApprox::LEguessApprox(const unsigned int _target_size,
			     Opponent _opp) :
    target_size(_target_size),
    opp(_opp),
    latest_entries(),
    min_unset(1)
{
    partial_lut[0] = 0;
    is_set[0] = true;
}


LEguessApprox::~LEguessApprox()
{
    partial_lut.clear();
    is_set.clear();
}


std::vector<IOpair> LEguessApprox::add_entry(const IOpair e)
{
    BinWord x = e.first, y = e.second;
    if ((x != 0) and (y == 0))
        opp.process_contradiction();
    else if ((x == 0) and (y != 0))
        opp.process_contradiction();
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
                        {
			     opp.process_contradiction();
			     break;
			}
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
        opp.process_contradiction();
    return std::vector<IOpair>();
}


bool LEguessApprox::is_entry_set(const BinWord x)
{
    return is_set[x];
}


unsigned int LEguessApprox::min_u()
{
    return min_unset;
}


bool LEguessApprox::complete()
{
    return (is_set.size() == target_size);
}


BinWord LEguessApprox::img(const BinWord x)
{
    if (is_set[x])
        return partial_lut[x];
    else
        throw std::runtime_error("Trying to access entry which doesn't exist");
}


Sbox LEguessApprox::lut()
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

LEguessIteratorApprox::LEguessIteratorApprox(const unsigned int _target_size,
					     const std::vector<IOpair> _constraints,
					     Opponent _opp) :
    target_size(_target_size),
    y(0),
    constraints(_constraints),
    opp(_opp),
    base_guess(target_size, opp),
    prepared_guess(target_size, opp)
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


bool LEguessIteratorApprox::prepare_successfully()
{
    y ++;
    if (y >= (target_size))
        return false;
    prepared_guess = LEguessApprox(base_guess);
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


LEguessApprox LEguessIteratorApprox::get_prepared_guess()
{
    return prepared_guess;
}



LEguessIteratorApprox LEguessIteratorApprox::deeper_guess_gen()
{
    std::vector<IOpair> deeper_constraints(constraints);
    deeper_constraints.push_back(IOpair(x, y));
    return LEguessIteratorApprox(target_size, deeper_constraints, opp);
}



// !SECTION! Actual algorithm LE

std::vector<Sbox> linear_equivalence_approx_cpp(const Sbox f,
						const Sbox g,
						const bool all_mappings,
						const unsigned int max_contradictions)
{
    check_length_cpp(f);
    check_length_cpp(g);
    std::vector<Sbox> result;
    // if f(0) == 0 then it is necessary that g(0) == 0...
    if ((f[0] == 0) and (g[0] != 0))
        return result;
    // ... and vice-versa
    else if ((f[0] != 0) and (g[0] == 0))
        return result;
    // We look for A,B such that f = B o g o A
    Sbox f_inv = inverse_cpp(f), g_inv = inverse_cpp(g);
    std::vector<IOpair> constraints;
    if (f[0] != 0)
    {                           // in this case, we know that B(g[0]) = f[0]
        constraints.push_back(IOpair(g[0], f[0]));
    }
    std::vector<LEguessIteratorApprox> b_generators;
    Opponent opp(max_contradictions);
    b_generators.push_back(LEguessIteratorApprox(f.size(), constraints, opp));
    while (true)
    {
        while ((b_generators.size() > 0) and (not b_generators.back().prepare_successfully()))
            b_generators.pop_back();
        if (b_generators.size() == 0)
                return result; // nothing more to be found, we have
                               // exhausted all possible guesses
        LEguessApprox
            a(f.size(), opp),
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
	opp.reset();
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
                Sbox A(f.size(), 0), B(b.lut());
                if (is_permutation_cpp(B))
                {
                    Sbox B_inv(inverse_cpp(B));
                    // A is rebuilt from scratch using that f = B o g o A
                    for (unsigned int x=0; x<B.size(); x++)
                        A[x] = g_inv[B_inv[f[x]]];
                    result.push_back(A);
                    result.push_back(B);
                    if (not all_mappings)
                        return result;
                }
                else
                    // inconclusive: we will guess one more entry
                    b_generators.push_back(b_generators.back().deeper_guess_gen());
                    
            }
        }
    }
}
