#include "./linear_representative.hpp"
#include <array>






// !SECTION! Lukas Stennes' fast linear representative, templated and generalized


typedef uint8_t u8;
typedef uint16_t u16;
typedef uint32_t u32;
typedef uint64_t u64;

template<typename set_t, typename int_type> struct tstate_t
{
        std::vector<int_type> A;
        std::vector<int_type> B;
        std::vector<int_type> R_S;
        set_t D_A;
        set_t D_B;
        set_t C_A;
        set_t C_B;
        set_t N_A;
        set_t N_B;
        set_t U_A;
        set_t U_B;
};


/*
 * set_t is the bit container
 * int_type is a type that can contain the sbox inputs
 */
template <typename set_t, typename int_type> class Set{
public:

    static void print(const std::string& name, const set_t& a);
    static inline int_type least_element(const set_t&a);
    static inline set_t intersect(const set_t&a, const set_t&b);
    static inline set_t unite(const set_t&a, const set_t&b);
    static inline u64 size(const set_t&a);
    static inline bool is_empty(const set_t&a);

    static inline set_t add_element(const set_t& a, const int_type elm);
    static inline set_t shift(const set_t& b, const int_type shift);
    static inline std::vector<int_type> get_elements(const set_t& a);

    static inline set_t init_empty(const u32 length);
    static inline set_t init_full(const u32 length);
    static inline set_t invert(const set_t& a, const u32 length);

    static inline bool elm_is_in_set(const int_type e, const set_t& a);

    static inline set_t shift_union(const set_t& a, const int_type shift)
    {
        set_t b = shift(a, shift);
        return unite(a, b);
    };

    static inline set_t setminus(const set_t& a, const set_t& b, const u32 length)
    {
        const set_t b_not = invert(b, length);
        return intersect(a, b_not);
    };
};

inline uint64_t permute_mask(uint64_t x, uint64_t mask, uint64_t len)
{
     return ((x & mask) << len) ^ ((x >> len) & mask);
}

template<typename int_type> class Set<std::vector<u64>, int_type>
{
    using set_t = std::vector<u64>;
public:
    static void print(const std::string& name, const set_t& a)
    {
        std::cout << name << ": 0b";
        for (int i = a.size() -1; i >= 0; i--)
        {
            for (int j = 63; j >= 0; j--)
            {
                u32 bit = (a[i] >> j) & 1;
                std::cout << bit;
            }
            std::cout << " ";
        }
        std::cout << std::endl;
    }

    static inline int_type least_element(const set_t& a)
    {
        for (u32 i = 0; i < a.size(); i++)
        {
            u64 current_chunk = a[i];
            if (current_chunk != 0)
            {
                u64 idx = __builtin_ctzll(current_chunk) + i * 64;
                return idx;
            }
        }

        //set is empty -- caller's fault
        std::cout << "CALLED LEAST ELEMENT ON EMPTY SET!" << std::endl;
        return 0;
    }

    static inline set_t intersect(const set_t& a, const set_t& b)
    {
        set_t result;
        result.reserve(a.size());
        for(u32 i = 0; i < a.size(); i++)
        {
            result.push_back(a[i] & b[i]);
        }
        return result;
    }

    static inline set_t unite(const set_t& a, const set_t& b)
    {
        set_t result;
        result.reserve(a.size());
        for(u64 i = 0; i < a.size(); i++)
        {
            result.push_back(a[i] | b[i]);
        }
        return result;
    }

    static inline u64 size(const set_t& a)
    {
        u64 count = 0;
        for(auto x : a)
        {
            count += __builtin_popcountll(x);
        }
        return count;
    }

    static inline bool is_empty(const set_t& a)
    {
        for(auto x : a)
        {
            if(x!=0)
                return false;
        }
        return true;
    }

    static inline set_t add_element(const set_t& a, const u64 elm)
    {
        //compute union of a and {elm}
        u64 index = elm / 64;
        set_t result = a;
        result[index] |= (1ul << (elm %64));
        return result;
    }

    static inline set_t shift(const set_t& b, const u64 shift)
    {
        set_t a = b;

        u64 high = shift >> 6;
        u64 msb = 1ul << cpp_msb(high);

        for( u64 bit = msb ; bit > 0; bit/=2){
            if (bit & high)
            {
            for(u64 i = 0; i < b.size(); i++)
                if( (i&bit) == 0)
                    std::swap(a[i],a[i+bit]);
            }
        }
        //compute a \oplus shift
        if ((shift >> 5) & 0x1)
        {
            for(int i = 0; i < b.size(); i++)
                a[i] = permute_mask(a[i], 0x00000000FFFFFFFF,32);
        }
        if ((shift >> 4) & 0x1)
        {
            for(int i = 0; i < b.size(); i++)
                a[i] = permute_mask(a[i], 0x0000FFFF0000FFFF,16);
        }
        if ((shift >> 3) & 0x1)
        {
            for(int i = 0; i < b.size(); i++)
                a[i] = permute_mask(a[i], 0x00FF00FF00FF00FF,8);
        }
        if ((shift >> 2) & 0x1)
        {
            for(u32 i = 0; i < b.size(); i++)
                a[i] = permute_mask(a[i], 0x0F0F0F0F0F0F0F0F,4);
        }
        if ((shift >> 1) & 0x1)
        {
            for(u32 i = 0; i < b.size(); i++)
                a[i] = permute_mask(a[i], 0x3333333333333333,2);
        }
        if (shift & 0x1)
        {
            for(u32 i = 0; i < b.size(); i++)
                a[i] = permute_mask(a[i], 0x5555555555555555,1);
        }
        return a;
    }

    static inline std::vector<int_type> get_elements(const set_t& a)
    {

        std::vector<u64> e;
        for (u64 i = 0; i < a.size(); i++)
        {
            u64 current_chunk = a[i];
            while (current_chunk != 0)
            {
                u64 idx = __builtin_ctzll(current_chunk) + i * 64;
                e.push_back(idx);
                current_chunk &= (current_chunk - 1);
            }
        }
        return e;
    }

    static inline set_t init_empty(const u32 len)
    {
        return set_t(len >> 6,0);
    }

    static inline set_t init_full(const u32 len)
    {
        u64 mask = -1ul;
        return set_t(len >> 6,mask);
    }

    static inline set_t invert(const set_t& a, const u32 len)
    {
        const set_t b = init_full(len);
        set_t result;
        result.reserve(a.size());
        for(u64 i = 0; i < a.size(); i++)
        {
            result.push_back(a[i] ^ b[i]);
        }
        return result;
    }

    static inline bool elm_is_in_set(const u64 e, const set_t& a)
    {
        u64 index = e >> 6;
        u64 field = 1ul << (e & 0x3F);
        return (a[index] & field) != 0;
    }

    static inline set_t shift_union(const set_t& a, const int_type s)
    {
        set_t b = shift(a, s);
        return unite(a, b);
    };

    static inline set_t setminus(const set_t& a, const set_t& b, const u32 length)
    {
        const set_t b_not = invert(b, length);
        return intersect(a, b_not);
    };


};

template <typename int_type, std::size_t s> class Set<std::array<u64, s>, int_type>{
    using set_t = std::array<u64, s>;
public:

    static void print(const std::string& name, const set_t& a)
    {
        std::cout << name << ": 0b";
        for (int i = s-1; i >=0; i--)
        {
            for (int j = 63; j >= 0; j--)
            {
                u32 bit = (a[i] >> j) & 1;
                std::cout << bit;
            }
            std::cout << " ";
        }
        std::cout << std::endl;
    }

    static inline int_type least_element(const set_t& a)
    {
        for (u32 i = 0; i < s; i++)
        {
            u64 current_chunk = a[i];
            if (current_chunk != 0)
            {
                int_type idx = __builtin_ctzll(current_chunk) + i * 64;
                return idx;
            }
        }

        //set is empty -- caller's fault
        std::cout << "CALLED LEAST ELEMENT ON EMPTY SET!" << std::endl;
        return 0;
    }

    static inline set_t intersect(const set_t& a, const set_t& b)
    {
        set_t result;
        for(u64 i = 0; i < s; i++)
        {
            result[i] = a[i] & b[i];
        }
        return result;
    }

    static inline set_t unite(const set_t& a, const set_t& b)
    {
        set_t result;
        for(u64 i = 0; i < s; i++)
        {
            result[i] = (a[i] | b[i]);
        }
        return result;
    }

    static inline u64 size(const set_t& a)
    {
        u64 count = 0;
        for(u64 i = 0; i < s; i++)
        {
            count += __builtin_popcountll(a[i]);
        }
        return count;
    }

    static inline bool is_empty(const set_t& a)
    {
        for(u64 i = 0; i < s; i++)
        {
            if(a[i]!=0)
                return false;
        }
        return true;
    }

    static inline set_t add_element(const set_t& a, const u64 elm)
    {
        //compute union of a and {elm}
        u64 index = s > 1 ? elm / 64 : 0;
        set_t result = a;
        result[index] |= (1ul << (elm %64));
        return result;
    }

    static inline set_t shift(const set_t& b, const int_type shift)
    {
        set_t a = b;

        if(s > 1)
        {
            u64 high = shift >> 6;
            u64 msb = 1ul << cpp_msb(high);

            for( u64 bit = msb ; bit > 0; bit/=2){
                if (bit & high)
                {
                for(u64 i = 0; i < b.size(); i++)
                    if( (i&bit) == 0)
                        std::swap(a[i],a[i+bit]);
                }
            }
        }
        //compute a \oplus shift
        if ((shift >> 5) & 0x1)
        {
            for(int i = 0; i < b.size(); i++)
                a[i] = permute_mask(a[i], 0x00000000FFFFFFFF,32);
        }
        if ((shift >> 4) & 0x1)
        {
            for(int i = 0; i < b.size(); i++)
                a[i] = permute_mask(a[i], 0x0000FFFF0000FFFF,16);
        }
        if ((shift >> 3) & 0x1)
        {
            for(int i = 0; i < b.size(); i++)
                a[i] = permute_mask(a[i], 0x00FF00FF00FF00FF,8);
        }
        if ((shift >> 2) & 0x1)
        {
            for(u32 i = 0; i < b.size(); i++)
                a[i] = permute_mask(a[i], 0x0F0F0F0F0F0F0F0F,4);
        }
        if ((shift >> 1) & 0x1)
        {
            for(u32 i = 0; i < b.size(); i++)
                a[i] = permute_mask(a[i], 0x3333333333333333,2);
        }
        if (shift & 0x1)
        {
            for(u32 i = 0; i < b.size(); i++)
                a[i] = permute_mask(a[i], 0x5555555555555555,1);
        }
        return a;
    }

    static inline std::vector<int_type> get_elements(const set_t& a)
    {

        std::vector<int_type> e;
        for (u64 i = 0; i < s; i++)
        {
            u64 current_chunk = a[i];
            while (current_chunk != 0)
            {
                int_type idx = __builtin_ctzll(current_chunk) + i * 64;
                e.push_back(idx);
                current_chunk &= (current_chunk - 1);
            }
        }
        return e;
    }

    static inline set_t init_empty(const u32 len)
    {
        return set_t{};
    }

    static inline set_t init_full(const u32 len)
    {

        u64 mask;
        if(s == 1 && len < 64)
            mask = (1ul << len) -1;
        else
            mask = -1ul;
        set_t set;
        for (auto &x : set)
            x = mask;
        return set;
    }

    static inline set_t invert(const set_t& a, const u32 len)
    {
        const set_t b = init_full(len);
        set_t result;
        for(u64 i = 0; i < s; i++)
        {
            result[i] = (a[i] ^ b[i]);
        }
        return result;
    }

    static inline bool elm_is_in_set(const u64 e, const set_t& a)
    {
        u64 index = s > 1 ? e >> 6 : 0;
        u64 field = 1ul << (e & 0x3F);
        return (a[index] & field) != 0;
    }

    static inline set_t shift_union(const set_t& a, const int_type sh)
    {
        set_t b = shift(a, sh);
        return unite(a, b);
    };

    static inline set_t setminus(const set_t& a, const set_t& b, const u32 length)
    {
        const set_t b_not = invert(b, length);
        return intersect(a, b_not);
    };

};

#if defined(__AVX2__) or defined(__ARM_NEON)
#ifdef __AVX2__
#include <immintrin.h>

using fastset_t = __m256i;
#elif defined(__ARM_NEON)
#include <arm_neon.h>

using fastset_t = uint64x2x2_t;
#endif

template <typename int_type> class Set<fastset_t, int_type>
{
    using set_t = fastset_t;
public:

    static void print(const std::string& name, const set_t& a)
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
    };

    static inline int_type least_element(const set_t&a)
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
                int_type idx = __builtin_ctzll(current_chunk) + i * 64;
                return idx;
            }
        }

        // set is empty -- caller's fault
        std::cout << "CALLED LEAST ELEMENT ON EMPTY SET!" << std::endl;
        return 0;
    };

    static inline set_t intersect(const set_t&a, const set_t&b)
    {
#ifdef __AVX2__
        return _mm256_and_si256(a, b);
#elif defined(__ARM_NEON)
        return {vandq_u64(a.val[0], b.val[0]), vandq_u64(a.val[1], b.val[1])};
#endif
    };

    static inline set_t unite(const set_t&a, const set_t&b)
    {
#ifdef __AVX2__
        return _mm256_or_si256(a, b);
#elif defined(__ARM_NEON)
        return {vorrq_u64(a.val[0], b.val[0]), vorrq_u64(a.val[1], b.val[1])};
#endif
    };

    static inline u64 size(const set_t&a)
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
    };

    static inline bool is_empty(const set_t&a)
    {
#ifdef __AVX2__
        return _mm256_testz_si256(a, a);
#elif defined(__ARM_NEON)
        auto tmp = vandq_u64(vceqzq_u64(a.val[0]), vceqzq_u64(a.val[1]));
        return (tmp[0] & tmp[1]) & 1;
#endif
    };

    static inline set_t add_element(const set_t& a, const int_type elm)
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
    };

    static inline set_t shift(const set_t& b, const int_type shift)
    {
        set_t a = b;
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
    };
    static inline std::vector<int_type> get_elements(const set_t& a)
    {
        std::vector<int_type> e;
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
                int_type idx = __builtin_ctzll(current_chunk) + i * 64;
                e.push_back(idx);
                current_chunk &= (current_chunk - 1);
            }
        }
        return e;
    };

    static inline set_t init_empty(const u32 length)
    {
#ifdef __AVX2__
        return _mm256_setzero_si256();
#elif defined(__ARM_NEON)
        return {vdupq_n_u64(0), vdupq_n_u64(0)};
#endif
    };
    static inline set_t init_full(const u32 len)
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
    };

    static inline set_t invert(const set_t& a, const u32 len)
    {
        const set_t b = init_full(len);
#ifdef __AVX2__
        return _mm256_xor_si256(a, b);
#elif defined(__ARM_NEON)
        return {veorq_u64(a.val[0], b.val[0]), veorq_u64(a.val[1], b.val[1])};
#endif
    };

    static inline bool elm_is_in_set(const int_type e, const set_t& a)
    {
        set_t b = init_empty(256);
        b            = add_element(b, e);
        b            = intersect(a, b);
        return !is_empty(b);
    };

    static inline set_t shift_union(const set_t& a, const int_type shift)
    {
        set_t b = shift(a, shift);
        return unite(a, b);
    };

    static inline set_t setminus(const set_t& a, const set_t& b, const u32 length)
    {
        const set_t b_not = invert(b, length);
        return intersect(a, b_not);
    };
};
#endif

template<typename int_type> bool update_linear(std::vector<int_type>& A, int_type new_x, const u32 length)
{
    int_type new_y = A[new_x];
    for (u32 i = 1; i < length; i++)
    {
        int_type e = A[i];
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

    // lexicographically compare R_S and R_S_best
template<typename int_type> bool is_greater(const std::vector<int_type>& R_S, const std::vector<int_type>& R_S_best, const u32 length)
{
    if ((R_S_best[0] == 0) && (R_S_best[1] == 0))
        return false;

    for (u32 x = 0; x < length; x++)
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

template<typename set_t, typename int_type> bool subroutine(const std::vector<int_type>& S, const std::vector<int_type>& S_inv, const tstate_t<set_t, int_type>& state, std::vector<int_type>& R_S_best, std::vector<int_type>& A_best, std::vector<int_type>& B_best, const u32 length)
{
    std::vector<int_type> A(state.A);
    std::vector<int_type> B(state.B);
    std::vector<int_type> R_S(state.R_S);

    set_t D_A = (state.D_A);
    set_t D_B = (state.D_B);
    set_t C_A = (state.C_A);
    set_t C_B = (state.C_B);
    set_t N_A = (state.N_A);
    set_t N_B = (state.N_B);
    set_t U_A = (state.U_A);
    set_t U_B = (state.U_B);

    while (!Set<set_t, int_type>::is_empty(N_A))
    {
        int_type x = Set<set_t, int_type>::least_element(N_A);
        int_type y = Set<set_t, int_type>::least_element(U_B);
        B[y] = S[A[x]];
        if (!update_linear<int_type>(B, y, length))
            return false;
        set_t D_B_new = Set<set_t, int_type>::shift(D_B, y);
        D_B           = Set<set_t, int_type>::unite(D_B, D_B_new);
        U_B           = Set<set_t, int_type>::setminus(U_B, D_B_new, length);
        set_t SoA_N_A = Set<set_t, int_type>::init_empty(length);
        for (int_type x : Set<set_t, int_type>::get_elements(N_A))
        {
            SoA_N_A = Set<set_t, int_type>::add_element(SoA_N_A, S[A[x]]);
        }
        set_t B_D_B_new = Set<set_t, int_type>::init_empty(length);
        for (int_type d : Set<set_t, int_type>::get_elements(D_B_new))
        {
            B_D_B_new = Set<set_t, int_type>::add_element(B_D_B_new, B[d]);
            if (Set<set_t, int_type>::elm_is_in_set(B[d], SoA_N_A))
            {
                C_B = Set<set_t, int_type>::add_element(C_B, d);
            }
            else
            {
                N_B = Set<set_t, int_type>::add_element(N_B, d);
            }
        }
        set_t C_A_new = Set<set_t, int_type>::init_empty(length);
        for (int_type x : Set<set_t, int_type>::get_elements(N_A))
        {
            if (Set<set_t, int_type>::elm_is_in_set(S[A[x]], B_D_B_new))
            {
                C_A_new = Set<set_t, int_type>::add_element(C_A_new, x);
            }
        }
        C_A = Set<set_t, int_type>::unite(C_A, C_A_new);
        N_A = Set<set_t, int_type>::setminus(N_A, C_A_new, length);
        for (int_type x : Set<set_t, int_type>::get_elements(C_A_new))
        {
            int_type y = 0;
            for (u32 i = 0; i < length; i++)
            {
                if (B[i] == S[A[x]])
                {
                    y = i;
                    break;
                }
            }
            R_S[x] = y;
        }
        if (is_greater<int_type>(R_S, R_S_best, length))
            return false;
        while (Set<set_t, int_type>::is_empty(N_A) && !Set<set_t, int_type>::is_empty(N_B))
        {
            int_type x = Set<set_t, int_type>::least_element(U_A);
            int_type y = Set<set_t, int_type>::least_element(N_B);
            A[x] = S_inv[B[y]];
            if (!update_linear<int_type>(A, x, length))
            {
                return false;
            }
            set_t D_A_new    = Set<set_t, int_type>::shift(D_A, x);
            D_A              = Set<set_t, int_type>::unite(D_A, D_A_new);
            U_A              = Set<set_t, int_type>::setminus(U_A, D_A_new, length);
            set_t SinvoB_N_B = Set<set_t, int_type>::init_empty(length);
            for (int_type y : Set<set_t, int_type>::get_elements(N_B))
            {
                SinvoB_N_B = Set<set_t, int_type>::add_element(SinvoB_N_B, S_inv[B[y]]);
            }
            set_t A_D_A_new = Set<set_t, int_type>::init_empty(length);
            for (int_type d : Set<set_t, int_type>::get_elements(D_A_new))
            {
                A_D_A_new = Set<set_t, int_type>::add_element(A_D_A_new, A[d]);
                if (Set<set_t, int_type>::elm_is_in_set(A[d], SinvoB_N_B))
                {
                    C_A = Set<set_t, int_type>::add_element(C_A, d);
                }
                else
                {
                    N_A = Set<set_t, int_type>::add_element(N_A, d);
                }
            }
            set_t C_B_new = Set<set_t, int_type>::init_empty(length);
            for (int_type y : Set<set_t, int_type>::get_elements(N_B))
            {
                if (Set<set_t, int_type>::elm_is_in_set(S_inv[B[y]], A_D_A_new))
                {
                    C_B_new = Set<set_t, int_type>::add_element(C_B_new, y);
                }
            }
            C_B = Set<set_t, int_type>::unite(C_B, C_B_new);
            N_B = Set<set_t, int_type>::setminus(N_B, C_B_new, length);
            for (int_type y : Set<set_t, int_type>::get_elements(C_B_new))
            {
                int_type x = 0;
                for (u32 i = 0; i < length; i++)
                {
                    if (A[i] == S_inv[B[y]])
                    {
                        x = i;
                        break;
                    }
                }
                R_S[x] = y;
            }
            if (is_greater<int_type>(R_S, R_S_best, length))
                return false;
        }
    }
    if (Set<set_t, int_type>::is_empty(U_A) && Set<set_t, int_type>::is_empty(U_B))
    {

        for (u32 i = 0; i < length; i++)
        {
            // new best
            R_S_best[i] = R_S[i];
            A_best[i] = A[i];
            B_best[i] = B[i];
        }
        return true;
    }
    else
    {
        int_type x         = Set<set_t, int_type>::least_element(U_A);
        set_t D_A_new      = Set<set_t, int_type>::shift(D_A, x);
        U_A                = Set<set_t, int_type>::setminus(U_A, D_A_new, length);
        D_A                = Set<set_t, int_type>::unite(D_A, D_A_new);
        N_A                = Set<set_t, int_type>::unite(N_A, D_A_new);
        bool flag = false;
        set_t Y            = Set<set_t, int_type>::init_full(length);
        set_t A_set        = Set<set_t, int_type>::init_empty(length);
        for (u32 i = 0; i < length; i++)
        {
            A_set = Set<set_t, int_type>::add_element(A_set, A[i]);
        }
        Y = Set<set_t, int_type>::setminus(Y, A_set, length);
        for (int_type y : Set<set_t, int_type>::get_elements(Y))
        {
            std::vector<int_type> A_next_guess(length);

            for (u32 i = 0; i < length; i++)
            {
                A_next_guess[i] = A[i];
            }
            A_next_guess[x] = y;
            if (!update_linear<int_type>(A_next_guess, x, length))
                continue;
            tstate_t<set_t, int_type> state_next;
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


            if (subroutine<set_t, int_type>(S, S_inv, state_next, R_S_best, A_best, B_best, length))
            {
                flag = true;
            }
        }

        return flag;
    }
}


template<typename set_t, typename int_type> std::array<std::vector<int_type>,3> compute_linear_representative(const std::vector<int_type>& sbox)
{
    u32 length = sbox.size();

    // variable for current best candidate
    std::vector<int_type> R_S_best(length, 0);
    std::vector<int_type> A_best(length, 0);
    std::vector<int_type> B_best(length, 0);

    // invert sbox
    std::vector<int_type> S_inv(length, 0);
    for (u32 x = 0; x < length; x++)
    {
        int_type y     = sbox[x];
        S_inv[y] = x;
    }

    // init the state of the algorithm
    tstate_t<set_t, int_type> state;
    state.A   = std::vector<int_type>(length, 0);
    state.B   = std::vector<int_type>(length, 0);
    state.R_S = std::vector<int_type>(length, 0);

    state.D_A = Set<set_t, int_type>::add_element(Set<set_t, int_type>::init_empty(length), 0);
    state.D_B = Set<set_t, int_type>::add_element(Set<set_t, int_type>::init_empty(length), 0);

    state.C_A = Set<set_t, int_type>::init_empty(length);
    state.C_B = Set<set_t, int_type>::init_empty(length);

    state.N_A = Set<set_t, int_type>::add_element(Set<set_t, int_type>::init_empty(length), 0);
    state.N_B = Set<set_t, int_type>::add_element(Set<set_t, int_type>::init_empty(length), 0);

    state.U_A = Set<set_t, int_type>::setminus(Set<set_t, int_type>::init_full(length), state.D_A, length);
    state.U_B = Set<set_t, int_type>::setminus(Set<set_t, int_type>::init_full(length), state.D_A, length);

    // init in special case S[0] == 0
    if (sbox[0] == 0)
    {
        state.C_A = Set<set_t, int_type>::add_element(Set<set_t, int_type>::init_empty(length), 0);
        state.C_B = Set<set_t, int_type>::add_element(Set<set_t, int_type>::init_empty(length), 0);

        state.N_A = Set<set_t, int_type>::init_empty(length);
        state.N_B = Set<set_t, int_type>::init_empty(length);
    }

    // compute linear representative recursively
    subroutine<set_t, int_type>(sbox, S_inv, state, R_S_best, A_best, B_best, length);

    return { R_S_best, A_best, B_best };
}


// !SECTION! The main function


/*
 * Current templates allow to use the same algorithm using various bit containers:
 * vectors, arrays, and 256-bit registers leveraging vectorial CPU instructions
 * We could have arrays of 256-bit registers, or even 512-bit registers or longer,
 * but given the limited gains compared to arrays of u64 this is probably not worth
 * the effort. Nevertheless, fixed-size arrays where benched ~2 times faster than
 * vectors, which explains the current approach.
 */

cpp_S_box cpp_le_class_representative(
    const cpp_S_box f,
    cpp_F2AffineMap & A,
    cpp_F2AffineMap & B)
{
    if(f.size() <= 256){
        std::vector<u8> ff;
        for(unsigned int i = 0; i < f.size(); i++){
            ff.push_back((u8) f[i]);
        }
        std::array<std::vector<u8>,3> result;
        if(f.size() <= 64)
            result = compute_linear_representative<std::array<u64,1>,u8>(ff);
        else if(f.size() == 128)
             result = compute_linear_representative<std::array<u64,2>,u8>(ff); // Could be faster with AVX/NEON/u128
        else
        {
#if defined(__AVX2__) or defined(__ARM_NEON)
             result = compute_linear_representative<fastset_t,u8>(ff);
#else
             result = compute_linear_representative<std::array<u64,4>,u8>(ff);
#endif
        }
        Lut linear_representative(f.size(), 0);
        Lut lut_A(f.size(), 0);
        Lut lut_B(f.size(), 0);

        for(unsigned int i = 0; i < f.size(); i++){
            linear_representative[i] = (BinWord)result[0][i];
            lut_A[i] = (BinWord)result[1][i];
            lut_B[i] = (BinWord)result[2][i];
        }
        A = cpp_F2AffineMap_from_lut(lut_A);
        B = cpp_F2AffineMap_from_lut(lut_B);
        return cpp_S_box(linear_representative);
    }
    std::vector<u64> ff;
    for(unsigned int i = 0; i < f.size(); i++){
        ff.push_back(f[i]);
    }
    std::array<std::vector<BinWord>,3> result;
    if(f.size() == 512)
         result = compute_linear_representative<std::array<u64,8>,u64>(ff); // Could be faster with AVX512
    else if(f.size() == 1024)
         result = compute_linear_representative<std::array<u64,16>,u64>(ff);
    else if(f.size() == 2048)
         result= compute_linear_representative<std::array<u64,32>,u64>(ff);
    else if(f.size() == 4096)
         result = compute_linear_representative<std::array<u64,64>,u64>(ff);
    else
        // Currently Set<std::vector<u64>, u64> does not support < 6 bit sboxes, due to simple init_full()
         result = compute_linear_representative<std::vector<u64>, u64>(ff);
    A = cpp_F2AffineMap_from_lut(result[1]);
    B = cpp_F2AffineMap_from_lut(result[2]);
    return cpp_S_box(result[0]);
}

