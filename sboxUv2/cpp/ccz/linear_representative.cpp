#include "./linear_representative.hpp"


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


Lut LEguess::lut()
{
    if (complete())
    {
        Lut result(target_size, 0);
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


// std::vector<Lut> linear_equivalence_cpp(const Lut f, const Lut g, bool all_mappings)
// {
//     check_length_cpp(f);
//     check_length_cpp(g);
//     std::vector<Lut> result;
//     // if f(0) == 0 then it is necessary that g(0) == 0...
//     if ((f[0] == 0) and (g[0] != 0))
//         return result;
//     // ... and vice-versa
//     else if ((f[0] != 0) and (g[0] == 0))
//         return result;
//     // We look for A,B such that f = B o g o A
//     Lut f_inv = inverse_cpp(f), g_inv = inverse_cpp(g);
//     std::vector<IOpair> constraints;
//     if (f[0] != 0)
//     {                           // in this case, we know that B(g[0]) = f[0]
//         constraints.push_back(IOpair(g[0], f[0]));
//     }
//     std::vector<LEguessIterator> b_generators;
//     b_generators.push_back(LEguessIterator(f.size(), constraints));
//     while (true)
//     {
//         while ((b_generators.size() > 0) and (not b_generators.back().prepare_successfully()))
//             b_generators.pop_back();
//         if (b_generators.size() == 0)
//                 return result; // nothing more to be found, we have
//                                // exhausted all possible guesses
//         LEguess
//             a(f.size()),
//             b = b_generators.back().get_prepared_guess();
        
//         std::vector<IOpair> just_added, to_treat, next_to_treat;
//         for (unsigned int x=0; x<f.size(); x++)
//             if (b.is_entry_set(x))
//                 to_treat.push_back(IOpair(x, b.img(x)));
//         bool
//             conclusion_reached = false, // true if and only if we have
//                                         // propagated the intial guess
//                                         // as far as possible.
//             contradiction_found = false; // true if and only if the
//                                          // initial guess leads to
//                                          // some inconsistency.
//         // propagating the guess on B
//         while (not conclusion_reached)
//         {
//             conclusion_reached = true;
//             try
//             {
//                 // B(y) = b_y => A(f_inv(b_y)) = g_inv(y)
//                 for (auto & entry : to_treat)
//                 {
//                     BinWord y = entry.first, b_y = entry.second;
//                     just_added = a.add_entry(IOpair(f_inv[b_y], g_inv[y]));
//                     next_to_treat.insert(next_to_treat.end(),
//                                          just_added.begin(),
//                                          just_added.end());
//                 }
//                 to_treat = next_to_treat;
//                 next_to_treat.clear();
//                 if (to_treat.size() > 0)
//                     conclusion_reached = false;
//                 // A(x) = a_x => B(g(a_x)) = f(x)
//                 for (auto & entry : to_treat)
//                 {
//                     BinWord x = entry.first, a_x = entry.second;
//                     just_added = b.add_entry(IOpair(g[a_x], f[x]));
//                     next_to_treat.insert(next_to_treat.end(),
//                                          just_added.begin(),
//                                          just_added.end());
//                 }
//                 to_treat = next_to_treat;
//                 next_to_treat.clear();
//                 if (to_treat.size() > 0)
//                     conclusion_reached = false;
//             }
//             catch (ContradictionFound & e)
//             {
//                 conclusion_reached = true;
//                 contradiction_found = true;
//             }
//         }
//         if (not contradiction_found)
//         {
//             if (not b.complete())
//                 // inconclusive: we will guess one more entry
//                 b_generators.push_back(b_generators.back().deeper_guess_gen());
//             else
//             {                   // we have everything we need
//                 Lut A(f.size(), 0), B(b.lut());
//                 if (is_permutation_cpp(B))
//                 {
//                     Lut B_inv(inverse_cpp(B));
//                     // A is rebuilt from scratch using that f = B o g o A
//                     for (unsigned int x=0; x<B.size(); x++)
//                         A[x] = g_inv[B_inv[f[x]]];
//                     result.push_back(A);
//                     result.push_back(B);
//                     if (not all_mappings)
//                         return result;
//                 }
//                 else
//                     // inconclusive: we will guess one more entry
//                     b_generators.push_back(b_generators.back().deeper_guess_gen());
                    
//             }
//         }
//     }
// }


// !SECTION! Finding the representative


// !SUBSECTION! The  unoptimized approach

LErepr::LErepr(const cpp_S_box _f) :
    f(_f),
    f_inv(_f.inverse()),
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


Lut LErepr::lut()
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



// !SUBSECTION! Lukas Stennes' fast linear representative 



typedef uint8_t u8;
typedef uint16_t u16;
typedef uint32_t u32;
typedef uint64_t u64;
typedef __uint128_t u128;

#ifdef __AVX2__
#include <immintrin.h>

using smallset_t = __m256i;
#elif defined(__ARM_NEON)
#include <arm_neon.h>

using smallset_t = uint64x2x2_t;
#endif


namespace
{
    void smallset_print(const std::string& name, const smallset_t& a)
    {
        u64 elements[4];
#ifdef __AVX2__

        elements[0] = _mm256_extract_epi64(a, 3);
        elements[1] = _mm256_extract_epi64(a, 2);
        elements[2] = _mm256_extract_epi64(a, 1);
        elements[3] = _mm256_extract_epi64(a, 0);
#elif defined(__ARM_NEON)
        elements[0] = a.val[1][1];
        elements[1] = a.val[1][0];
        elements[2] = a.val[0][1];
        elements[3] = a.val[0][0];
#endif
        std::cout << name << ": 0b";
        for (u32 i = 0; i < 4; i++)
        {
            for (int j = 63; j >= 0; j--)
            {
                u32 bit = (elements[i] >> j) & 1;
                std::cout << bit;
            }
            std::cout << " ";
        }
        std::cout << std::endl;
    }

    u8 smallset_least_element(const smallset_t& a)
    {
        u64 chunks[4];
#ifdef __AVX2__
        chunks[0] = _mm256_extract_epi64(a, 0);
        chunks[1] = _mm256_extract_epi64(a, 1);
        chunks[2] = _mm256_extract_epi64(a, 2);
        chunks[3] = _mm256_extract_epi64(a, 3);
#elif defined(__ARM_NEON)
        chunks[0] = a.val[0][0];
        chunks[1] = a.val[0][1];
        chunks[2] = a.val[1][0];
        chunks[3] = a.val[1][1];
#endif
        for (u32 i = 0; i < 4; i++)
        {
            u64 current_chunk = chunks[i];
            if (current_chunk != 0)
            {
                u8 idx = __builtin_ctzll(current_chunk) + i * 64;
                return idx;
            }
        }

        // set is empty -- caller's fault
        std::cout << "CALLED LEAST ELEMENT ON EMPTY SET!" << std::endl;
        return 0;
    }

    inline smallset_t smallset_intersect(const smallset_t& a, const smallset_t& b)
    {
#ifdef __AVX2__
        return _mm256_and_si256(a, b);
#elif defined(__ARM_NEON)
        return {vandq_u64(a.val[0], b.val[0]), vandq_u64(a.val[1], b.val[1])};
#endif
    }

    inline smallset_t smallset_union(const smallset_t& a, const smallset_t& b)
    {
#ifdef __AVX2__
        return _mm256_or_si256(a, b);
#elif defined(__ARM_NEON)
        return {vorrq_u64(a.val[0], b.val[0]), vorrq_u64(a.val[1], b.val[1])};
#endif
    }

    inline u16 smallset_size(const smallset_t& a)
    {
        u16 count = 0;
#ifdef __AVX2__
        u64 chunk = _mm256_extract_epi64(a, 0);
        count += __builtin_popcountll(chunk);
        chunk = _mm256_extract_epi64(a, 1);
        count += __builtin_popcountll(chunk);
        chunk = _mm256_extract_epi64(a, 2);
        count += __builtin_popcountll(chunk);
        chunk = _mm256_extract_epi64(a, 3);
        count += __builtin_popcountll(chunk);
#elif defined(__ARM_NEON)
        count += __builtin_popcountll(a.val[0][0]);
        count += __builtin_popcountll(a.val[0][1]);
        count += __builtin_popcountll(a.val[1][0]);
        count += __builtin_popcountll(a.val[1][1]);
#endif
        return count;
    }

    inline bool smallset_is_empty(const smallset_t& a)
    {
#ifdef __AVX2__
        return _mm256_testz_si256(a, a);
#elif defined(__ARM_NEON)
        auto tmp = vandq_u64(vceqzq_u64(a.val[0]), vceqzq_u64(a.val[1]));
        return (tmp[0] & tmp[1]) & 1;
#endif
    }

    smallset_t smallset_add_element(const smallset_t& a, const u8 elm)
    {
        // compute union of a and {elm}
        u32 index = elm / 64;
#ifdef __AVX2__
        u64 mask[4]   = {0};
        mask[index]   = (u64)1 << (elm % 64);
        __m256i _mask = _mm256_set_epi64x(mask[3], mask[2], mask[1], mask[0]);
        return _mm256_or_si256(a, _mask);
#elif defined(__ARM_NEON)
        u64 mask[2]     = {0};
        mask[index & 1] = (u64)1 << (elm % 64);
        auto _mask      = vld1q_u64(mask);
        if (index < 2)
        {
            return {vorrq_u64(a.val[0], _mask), a.val[1]};
        }
        else
        {
            return {a.val[0], vorrq_u64(a.val[1], _mask)};
        }
#endif
    }

    smallset_t smallset_shift(const smallset_t& b, const u8 shift)
    {
        auto a = b;
        // compute a \oplus shift
        if ((shift >> 7) & 0x1)
        {
#ifdef __AVX2__
            a = _mm256_permute2x128_si256(a, a, 1);
#elif defined(__ARM_NEON)
            a.val[0] = b.val[1];
            a.val[1] = b.val[0];
#endif
        }
        if ((shift >> 6) & 0x1)
        {
#ifdef __AVX2__
            a = _mm256_permute4x64_epi64(a, _MM_SHUFFLE(2, 3, 0, 1));
#elif defined(__ARM_NEON)
            a.val[0] = vextq_u64(a.val[0], a.val[0], 1);
            a.val[1] = vextq_u64(a.val[1], a.val[1], 1);
#endif
        }
        if ((shift >> 5) & 0x1)
        {
#ifdef __AVX2__
            a = _mm256_shuffle_epi32(a, _MM_SHUFFLE(2, 3, 0, 1));
#elif defined(__ARM_NEON)
            a.val[0] = vrev64q_u32(a.val[0]);
            a.val[1] = vrev64q_u32(a.val[1]);

#endif
        }
        if ((shift >> 4) & 0x1)
        {
#ifdef __AVX2__
            a = _mm256_shufflelo_epi16(a, _MM_SHUFFLE(2, 3, 0, 1));
            a = _mm256_shufflehi_epi16(a, _MM_SHUFFLE(2, 3, 0, 1));
#elif defined(__ARM_NEON)
            a.val[0] = vrev64q_u16(a.val[0]);
            a.val[0] = vrev64q_u32(a.val[0]);
            a.val[1] = vrev64q_u16(a.val[1]);
            a.val[1] = vrev64q_u32(a.val[1]);
#endif
        }
        if ((shift >> 3) & 0x1)
        {
#ifdef __AVX2__
            const __m256i mask = _mm256_set_epi8(14, 15, 12, 13, 10, 11, 8, 9, 6, 7, 4, 5, 2, 3, 0, 1, 14, 15, 12, 13, 10, 11, 8, 9, 6, 7, 4, 5, 2, 3, 0, 1);
            a                  = _mm256_shuffle_epi8(a, mask);
#elif defined(__ARM_NEON)
            a.val[0] = vrev64q_u8(a.val[0]);
            a.val[0] = vrev64q_u16(a.val[0]);
            a.val[1] = vrev64q_u8(a.val[1]);
            a.val[1] = vrev64q_u16(a.val[1]);
#endif
        }
        if ((shift >> 2) & 0x1)
        {
#ifdef __AVX2__
            const __m256i mask_high = _mm256_set1_epi8((char)0xF0);
            const __m256i mask_low  = _mm256_set1_epi8(0x0F);
            const __m256i high      = _mm256_and_si256(a, mask_high);
            const __m256i low       = _mm256_and_si256(a, mask_low);
            a                       = _mm256_or_si256(_mm256_srli_epi16(high, 4), _mm256_slli_epi16(low, 4));
#elif defined(__ARM_NEON)
            const auto mask_high = vdupq_n_u64(0xF0F0F0F0F0F0F0F0);
            const auto mask_low  = vdupq_n_u64(0x0F0F0F0F0F0F0F0F);

            for (u32 i = 0; i < 2; i++)
            {
                const auto high = vandq_u64(a.val[i], mask_high);
                const auto low  = vandq_u64(a.val[i], mask_low);

                a.val[i] = vorrq_u64(vshrq_n_u64(high, 4), vshlq_n_u64(low, 4));
            }
#endif
        }
        if ((shift >> 1) & 0x1)
        {
#ifdef __AVX2__
            const __m256i mask_high = _mm256_set1_epi8((char)0xCC);
            const __m256i mask_low  = _mm256_set1_epi8(0x33);
            const __m256i high      = _mm256_and_si256(a, mask_high);
            const __m256i low       = _mm256_and_si256(a, mask_low);
            a                       = _mm256_or_si256(_mm256_srli_epi16(high, 2), _mm256_slli_epi16(low, 2));
#elif defined(__ARM_NEON)
            const auto mask_high = vdupq_n_u64(0xCCCCCCCCCCCCCCCC);
            const auto mask_low  = vdupq_n_u64(0x3333333333333333);

            for (u32 i = 0; i < 2; i++)
            {
                const auto high = vandq_u64(a.val[i], mask_high);
                const auto low  = vandq_u64(a.val[i], mask_low);

                a.val[i] = vorrq_u64(vshrq_n_u64(high, 2), vshlq_n_u64(low, 2));
            }
#endif
        }
        if (shift & 0x1)
        {
#ifdef __AVX2__
            const __m256i mask_high = _mm256_set1_epi8((char)0xAA);
            const __m256i mask_low  = _mm256_set1_epi8(0x55);
            const __m256i high      = _mm256_and_si256(a, mask_high);
            const __m256i low       = _mm256_and_si256(a, mask_low);
            a                       = _mm256_or_si256(_mm256_srli_epi16(high, 1), _mm256_slli_epi16(low, 1));
#elif defined(__ARM_NEON)
            const auto mask_high = vdupq_n_u64(0xAAAAAAAAAAAAAAAA);
            const auto mask_low  = vdupq_n_u64(0x5555555555555555);

            for (u32 i = 0; i < 2; i++)
            {
                const auto high = vandq_u64(a.val[i], mask_high);
                const auto low  = vandq_u64(a.val[i], mask_low);

                a.val[i] = vorrq_u64(vshrq_n_u64(high, 1), vshlq_n_u64(low, 1));
            }
#endif
        }
        return a;
    }

    smallset_t smallset_shift_union(const smallset_t& a, const u8 shift)
    {
        smallset_t b = smallset_shift(a, shift);
        return smallset_union(a, b);
    }

    std::vector<u8> smallset_get_elements(const smallset_t& a)
    {
        std::vector<u8> e;
        u64 chunks[4];
#ifdef __AVX2__
        chunks[0] = _mm256_extract_epi64(a, 0);
        chunks[1] = _mm256_extract_epi64(a, 1);
        chunks[2] = _mm256_extract_epi64(a, 2);
        chunks[3] = _mm256_extract_epi64(a, 3);
#elif defined(__ARM_NEON)
        chunks[0] = a.val[0][0];
        chunks[1] = a.val[0][1];
        chunks[2] = a.val[1][0];
        chunks[3] = a.val[1][1];
#endif
        for (u32 i = 0; i < 4; i++)
        {
            u64 current_chunk = chunks[i];
            while (current_chunk != 0)
            {
                u8 idx = __builtin_ctzll(current_chunk) + i * 64;
                e.push_back(idx);
                current_chunk &= (current_chunk - 1);
            }
        }
        return e;
    }

    inline smallset_t smallset_init_empty()
    {
#ifdef __AVX2__
        return _mm256_setzero_si256();
#elif defined(__ARM_NEON)
        return {vdupq_n_u64(0), vdupq_n_u64(0)};
#endif
    }

    inline smallset_t smallset_init_full(const u32 len)
    {
        // N must be in {256, 128, 64, 32, 16, 8}
        if (len == 256)
        {
#ifdef __AVX2__
            return _mm256_set_epi64x(0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF);
#elif defined(__ARM_NEON)
            return {vdupq_n_u64(0xFFFFFFFFFFFFFFFF), vdupq_n_u64(0xFFFFFFFFFFFFFFFF)};
#endif
        }
        else if (len == 128)
        {
#ifdef __AVX2__
            return _mm256_set_epi64x(0, 0, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF);
#elif defined(__ARM_NEON)
            return {vdupq_n_u64(0xFFFFFFFFFFFFFFFF), vdupq_n_u64(0)};
#endif
        }
        else if (len == 64)
        {
#ifdef __AVX2__
            return _mm256_set_epi64x(0, 0, 0, 0xFFFFFFFFFFFFFFFF);
#elif defined(__ARM_NEON)
            auto tmp = vdupq_n_u64(0);
            return {vsetq_lane_u64(0xFFFFFFFFFFFFFFFF, tmp, 0), vdupq_n_u64(0)};
#endif
        }
        else if (len == 32)
        {
#ifdef __AVX2__
            return _mm256_set_epi64x(0, 0, 0, 0xFFFFFFFF);
#elif defined(__ARM_NEON)
            auto tmp = vdupq_n_u64(0);
            return {vsetq_lane_u32(0xFFFFFFFF, tmp, 0), vdupq_n_u64(0)};
#endif
        }
        else if (len == 16)
        {
#ifdef __AVX2__
            return _mm256_set_epi64x(0, 0, 0, 0xFFFF);
#elif defined(__ARM_NEON)
            auto tmp = vdupq_n_u64(0);
            return {vsetq_lane_u16(0xFFFF, tmp, 0), vdupq_n_u64(0)};
#endif
        }
        else if (len == 8)
        {
#ifdef __AVX2__
            return _mm256_set_epi64x(0, 0, 0, 0xFF);
#elif defined(__ARM_NEON)
            auto tmp = vdupq_n_u64(0);
            return {vsetq_lane_u8(0xFF, tmp, 0), vdupq_n_u64(0)};
#endif
        }
        else
        {
#ifdef __AVX2__
            return _mm256_set_epi64x(0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF);
#elif defined(__ARM_NEON)
            return {vdupq_n_u64(0xFFFFFFFFFFFFFFFF), vdupq_n_u64(0xFFFFFFFFFFFFFFFF)};
#endif
        }
    }

    inline smallset_t smallset_invert(const smallset_t& a, const u32 len)
    {
        const smallset_t b = smallset_init_full(len);
#ifdef __AVX2__
        return _mm256_xor_si256(a, b);
#elif defined(__ARM_NEON)
        return {veorq_u64(a.val[0], b.val[0]), veorq_u64(a.val[1], b.val[1])};
#endif
    }

    inline smallset_t smallset_setminus(const smallset_t& a, const smallset_t& b, const u32 len)
    {
        const smallset_t b_not = smallset_invert(b, len);
        return smallset_intersect(a, b_not);
    }

    bool smallset_elm_is_in_set(const u8 e, const smallset_t& a)
    {
        smallset_t b = smallset_init_empty();
        b            = smallset_add_element(b, e);
        b            = smallset_intersect(a, b);
        return !smallset_is_empty(b);
    }
    // END OF SMALL SET //

    // state of the linear_representative algorithm
    typedef struct
    {
        std::vector<u8> A;
        std::vector<u8> B;
        std::vector<u8> R_S;
        smallset_t D_A;
        smallset_t D_B;
        smallset_t C_A;
        smallset_t C_B;
        smallset_t N_A;
        smallset_t N_B;
        smallset_t U_A;
        smallset_t U_B;
    } state_t;

    // lexicographically compare R_S and R_S_best
    bool is_greater(const std::vector<u8>& R_S, const std::vector<u8>& R_S_best, const u32 len)
    {
        if ((R_S_best[0] == 0) && (R_S_best[1] == 0))
            return false;

        for (u32 x = 0; x < len; x++)
        {
            // special case: R_S[x] not defined (=> 0) and R_S_best[x] = 0
            // works out with this
            if (R_S[x] > R_S_best[x])
                return true;
            if (R_S[x] < R_S_best[x])
                return false;
        }
        // can happen if there are self equivalences (?)
        return false;
    }

    bool update_linear(std::vector<u8>& A, u8 new_x, const u32 len)
    {
        u8 new_y = A[new_x];
        for (u32 i = 1; i < len; i++)
        {
            u8 e = A[i];
            if (e == 0)
                continue;
            else if (A[new_x ^ i] == 0)
                A[new_x ^ i] = e ^ new_y;
            else if (A[new_x ^ i] != (e ^ new_y))
            {
                return false;
            }
        }
        return true;
    }

    bool subroutine(const std::vector<u8>& S, const std::vector<u8>& S_inv, const state_t& state, std::vector<u8>& R_S_best, const u32 len)
    {
        std::vector<u8> A(state.A);
        std::vector<u8> B(state.B);
        std::vector<u8> R_S(state.R_S);

        smallset_t D_A = (state.D_A);
        smallset_t D_B = (state.D_B);
        smallset_t C_A = (state.C_A);
        smallset_t C_B = (state.C_B);
        smallset_t N_A = (state.N_A);
        smallset_t N_B = (state.N_B);
        smallset_t U_A = (state.U_A);
        smallset_t U_B = (state.U_B);

        while (!smallset_is_empty(N_A))
        {
            u8 x = smallset_least_element(N_A);
            u8 y = smallset_least_element(U_B);
            B[y] = S[A[x]];
            if (!update_linear(B, y, len))
                return false;
            smallset_t D_B_new = smallset_shift(D_B, y);
            D_B                = smallset_union(D_B, D_B_new);
            U_B                = smallset_setminus(U_B, D_B_new, len);
            smallset_t SoA_N_A = smallset_init_empty();
            for (u8 x : smallset_get_elements(N_A))
            {
                SoA_N_A = smallset_add_element(SoA_N_A, S[A[x]]);
            }
            smallset_t B_D_B_new = smallset_init_empty();
            for (u8 d : smallset_get_elements(D_B_new))
            {
                B_D_B_new = smallset_add_element(B_D_B_new, B[d]);
                if (smallset_elm_is_in_set(B[d], SoA_N_A))
                {
                    C_B = smallset_add_element(C_B, d);
                }
                else
                {
                    N_B = smallset_add_element(N_B, d);
                }
            }
            smallset_t C_A_new = smallset_init_empty();
            for (u8 x : smallset_get_elements(N_A))
            {
                if (smallset_elm_is_in_set(S[A[x]], B_D_B_new))
                {
                    C_A_new = smallset_add_element(C_A_new, x);
                }
            }
            C_A = smallset_union(C_A, C_A_new);
            N_A = smallset_setminus(N_A, C_A_new, len);
            for (u8 x : smallset_get_elements(C_A_new))
            {
                u8 y = 0;
                for (u32 i = 0; i < len; i++)
                {
                    if (B[i] == S[A[x]])
                    {
                        y = i;
                        break;
                    }
                }
                R_S[x] = y;
            }
            if (is_greater(R_S, R_S_best, len))
                return false;
            while (smallset_is_empty(N_A) && !smallset_is_empty(N_B))
            {
                u8 x = smallset_least_element(U_A);
                u8 y = smallset_least_element(N_B);
                A[x] = S_inv[B[y]];
                if (!update_linear(A, x, len))
                {
                    return false;
                }
                smallset_t D_A_new    = smallset_shift(D_A, x);
                D_A                   = smallset_union(D_A, D_A_new);
                U_A                   = smallset_setminus(U_A, D_A_new, len);
                smallset_t SinvoB_N_B = smallset_init_empty();
                for (u8 y : smallset_get_elements(N_B))
                {
                    SinvoB_N_B = smallset_add_element(SinvoB_N_B, S_inv[B[y]]);
                }
                smallset_t A_D_A_new = smallset_init_empty();
                for (u8 d : smallset_get_elements(D_A_new))
                {
                    A_D_A_new = smallset_add_element(A_D_A_new, A[d]);
                    if (smallset_elm_is_in_set(A[d], SinvoB_N_B))
                    {
                        C_A = smallset_add_element(C_A, d);
                    }
                    else
                    {
                        N_A = smallset_add_element(N_A, d);
                    }
                }
                smallset_t C_B_new = smallset_init_empty();
                for (u8 y : smallset_get_elements(N_B))
                {
                    if (smallset_elm_is_in_set(S_inv[B[y]], A_D_A_new))
                    {
                        C_B_new = smallset_add_element(C_B_new, y);
                    }
                }
                C_B = smallset_union(C_B, C_B_new);
                N_B = smallset_setminus(N_B, C_B_new, len);
                for (u8 y : smallset_get_elements(C_B_new))
                {
                    u8 x = 0;
                    for (u32 i = 0; i < len; i++)
                    {
                        if (A[i] == S_inv[B[y]])
                        {
                            x = i;
                            break;
                        }
                    }
                    R_S[x] = y;
                }
                if (is_greater(R_S, R_S_best, len))
                    return false;
            }
        }
        if (smallset_is_empty(U_A) && smallset_is_empty(U_B))
        {
            for (u32 i = 0; i < len; i++)
            {
                // new best
                R_S_best[i] = R_S[i];
            }
            return true;
        }
        else
        {
            u8 x               = smallset_least_element(U_A);
            smallset_t D_A_new = smallset_shift(D_A, x);
            U_A                = smallset_setminus(U_A, D_A_new, len);
            D_A                = smallset_union(D_A, D_A_new);
            N_A                = smallset_union(N_A, D_A_new);
            bool flag          = false;
            smallset_t Y       = smallset_init_full(len);
            smallset_t A_set   = smallset_init_empty();
            for (u32 i = 0; i < len; i++)
            {
                A_set = smallset_add_element(A_set, A[i]);
            }
            Y = smallset_setminus(Y, A_set, len);
            for (u8 y : smallset_get_elements(Y))
            {
                std::vector<u8> A_next_guess(len);

                for (u32 i = 0; i < len; i++)
                {
                    A_next_guess[i] = A[i];
                }
                A_next_guess[x] = y;
                if (!update_linear(A_next_guess, x, len))
                    continue;
                state_t state_next;
                state_next.A   = A_next_guess;
                state_next.B   = B;
                state_next.R_S = R_S;
                state_next.D_A = D_A;
                state_next.D_B = D_B;
                state_next.C_A = C_A;
                state_next.C_B = C_B;
                state_next.N_A = N_A;
                state_next.N_B = N_B;
                state_next.U_A = U_A;
                state_next.U_B = U_B;

                if (subroutine(S, S_inv, state_next, R_S_best, len))
                {
                    flag = true;
                }
            }

            return flag;
        }
    }
}    // namespace

std::vector<u8> compute_linear_representative(const std::vector<u8>& sbox)
{
    u32 len = sbox.size();

    // variable for current best candidate
    std::vector<u8> R_S_best(len, 0);

    // invert sbox
    std::vector<u8> S_inv(len, 0);
    for (u32 x = 0; x < len; x++)
    {
        u8 y     = sbox[x];
        S_inv[y] = x;
    }

    // init the state of the algorithm
    state_t state;
    state.A   = std::vector<u8>(len, 0);
    state.B   = std::vector<u8>(len, 0);
    state.R_S = std::vector<u8>(len, 0);

    state.D_A = smallset_add_element(smallset_init_empty(), 0);
    state.D_B = smallset_add_element(smallset_init_empty(), 0);

    state.C_A = smallset_init_empty();
    state.C_B = smallset_init_empty();

    state.N_A = smallset_add_element(smallset_init_empty(), 0);
    state.N_B = smallset_add_element(smallset_init_empty(), 0);

    state.U_A = smallset_setminus(smallset_init_full(len), state.D_A, len);
    state.U_B = smallset_setminus(smallset_init_full(len), state.D_A, len);

    // init in special case S[0] == 0
    if (sbox[0] == 0)
    {
        state.C_A = smallset_add_element(smallset_init_empty(), 0);
        state.C_B = smallset_add_element(smallset_init_empty(), 0);

        state.N_A = smallset_init_empty();
        state.N_B = smallset_init_empty();
    }

    // compute linear representative recursively
    subroutine(sbox, S_inv, state, R_S_best, len);

    return R_S_best;
}


// !SUBSECTION! The main function

cpp_S_box cpp_le_class_representative(const cpp_S_box f, const uint64_t fast)
{
    if(f.size() <= 256 and fast){
      std::vector<u8> ff;
      for(unsigned int i = 0; i < f.size(); i++){
        ff.push_back((u8) f[i]);
      }
      std::vector<u8> lr = compute_linear_representative(ff);
      Lut linear_representative(f.size(), 0);
      for(unsigned int i = 0; i < f.size(); i++){
          linear_representative[i] = (BinWord)lr[i];
      }
      return linear_representative;
    }
    LErepr repr(f);
    repr.initialize();
    return cpp_S_box(repr.lut());
}
