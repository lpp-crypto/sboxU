#ifndef _BIG_VECTORS_
#define _BIG_VECTORS_

#define MSB_OF_ZERO -1

#define BoolBlock uint64_t
#define BLOCK_SIZE (sizeof(BoolBlock) * 8)
#define BLOCK_INDEX(x) ((x) >> 6)
#define BLOCK_POS(x) ((x) & 0x3F)


#include <valarray>

// !SECTION! Big vectors themselves


class cpp_BigF2Vector
{
private:
    int msb;
    unsigned int total_length;
    
public:
    std::valarray<BoolBlock> content;

    
    // !SUBSECTION! constructors 

    cpp_BigF2Vector() :
        content(),
        msb(MSB_OF_ZERO),
        total_length(0)
    {}

        
    cpp_BigF2Vector(unsigned int _total_length) :
        content(),
        msb(MSB_OF_ZERO),
        total_length(_total_length)
    {
        if ((total_length % BLOCK_SIZE) == 0)
            content = std::valarray<BoolBlock>(total_length/BLOCK_SIZE);
        else
            content = std::valarray<BoolBlock>(total_length/BLOCK_SIZE+1);
    }

    
    cpp_BigF2Vector(
        std::valarray<BoolBlock> _content,
        unsigned int _total_length
        ) :
        content(_content),
        msb(MSB_OF_ZERO),
        total_length(_total_length)
    {
        set_msb();
    }

    
    void set_msb()
    {
        for(unsigned int i=content.size()-1; i>0; i--)
            if (content[i])
            {
                msb = i*BLOCK_SIZE + cpp_msb(content[i]);
                return;
            }
        if (content[0])
            msb = cpp_msb(content[0]);
        else
            msb = MSB_OF_ZERO;
    }

    
    // !SUBSECTION! Reading the state


    inline unsigned int size() const
    {
        return total_length;
    }
    
    
    inline unsigned int n_blocks() const
    {
        return content.size();
    }

    /** Returns 1 if the bit at position index is set to 1, 0 otherwise.

        @param index The position at which the bit must be set to 1.

        @return true if and only if the bit at position index is set to 1.
     */
    inline bool is_set(const unsigned int index) const
    {
        return (content[BLOCK_INDEX(index)] >> BLOCK_POS(index)) & 1;
    }
    
    inline unsigned int get_msb() const
    {
        return msb;
    }
      

    inline bool is_zero() const
    {
        return (msb == MSB_OF_ZERO);
    }


    inline bool is_non_zero() const
    {
        return (msb != MSB_OF_ZERO);
    }


    std::string to_string() const
    {
        std::stringstream result;
        result << std::setw(2*sizeof(BoolBlock)) << std::hex ;
        for(auto x_i : content)
            result << x_i << "  ";
        return result.str();
    }
    

    Bytearray to_bytes() const
    {
        Bytearray result;
        result.reserve(sizeof(BoolBlock) * content.size());
        for(auto v : content)
            for(unsigned int i=0; i<BLOCK_SIZE; i+=8)
                result.push_back((v >> i) & 0xFF);
        return result;
    }

    
    Bytearray to_bits() const
    {
        Bytearray result;
        result.reserve(8 * sizeof(BoolBlock) * content.size());
        for(auto v : content)
            for(unsigned int i=0; i<BLOCK_SIZE; i++)
                result.push_back((v >> i) & 1);
        return result;
    }

    
    // !SUBSECTION! modifying the state

    /** Ensures that the bit at position index is set to 1.

        Like is_set, its implementation relies on bit-fiddling.

        @param index The position at which the bit must be set to 1.
     */
    inline void set_to_1(const unsigned int index)
    {
        if (msb == MSB_OF_ZERO)
        {
            content[BLOCK_INDEX(index)] = ((BoolBlock)1 << (BLOCK_POS(index)));
            msb = index;
        }
        else
        {
            content[BLOCK_INDEX(index)] |= ((BoolBlock)1 << (BLOCK_POS(index)));
            if (index > msb)
                msb = index;
        }
    }

    inline void operator^=(const cpp_BigF2Vector & x)
    {
        content ^= x.content;
        set_msb();
    }
};


// !SECTION! Operators and helper functions

inline cpp_BigF2Vector operator^(const cpp_BigF2Vector & x,
                                 const cpp_BigF2Vector & y)
{
    return cpp_BigF2Vector(x.content ^ y.content, x.size());
}


inline bool operator==(const cpp_BigF2Vector & x,
                       const cpp_BigF2Vector & y)
{
    if (x.get_msb() == y.get_msb())
        if (x.is_zero())
            return true;
        else
        {
            for(unsigned int i=0; i<=BLOCK_INDEX(x.get_msb()); i++)
                if (x.content[i] ^ y.content[i])
                    return false;
            return true;
        }
    else
        return false;
}


inline bool operator<(const cpp_BigF2Vector & x,
                      const cpp_BigF2Vector & y)
{
    if (x.get_msb() != y.get_msb())
        return (x.get_msb() < y.get_msb());
    else
    {
        for (int i=x.n_blocks()-1; i >= 0; i--)
        {
            if (x.content[i] < y.content[i])
                return true;
            else if (x.content[i] > y.content[i])
                return false;
        }
        // if we reach this point, then they are equal. Thus, the
        // strict comparison is false.
        return false;
    }
}

#endif
