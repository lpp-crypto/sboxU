#ifndef _BIG_VECTORS_
#define _BIG_VECTORS_

#define MSB_OF_ZERO 0xffffffff // 2**32-1
#define BLOCK_SIZE (sizeof(BinWord)*8)


// !SECTION! Big vectors themselves


class cpp_BigF2Vector
{
private:
    unsigned int msb;
    unsigned int total_length;
    
public:
    std::vector<BinWord> content;

    
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
        content.assign(1+total_length/BLOCK_SIZE, 0);
    }

    
    cpp_BigF2Vector(
        std::vector<BinWord> _content,
        unsigned int _total_length
        ) :
        content(_content.begin(), _content.end()),
        msb(MSB_OF_ZERO),
        total_length(_total_length)
    {
        if (content.size() != (1+ total_length/BLOCK_SIZE))
            throw std::runtime_error("mismatched length between inputs");
        set_msb();
    }

    
    void set_msb()
    {
        unsigned int cursor = 0;
        for (unsigned int i=0; i<content.size(); i++)
        {
            BinWord x_i = content[i];
            if (x_i != 0)
                msb = cursor + cpp_msb(x_i);
            cursor += BLOCK_SIZE;
        }
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

    
    inline bool is_set(const unsigned int index) const
    {
        unsigned int
            cursor = index / BLOCK_SIZE,
            pos    = index % BLOCK_SIZE;

        // std::cout << std::dec << index << std::endl;
        // std::cout << std::dec << BLOCK_SIZE << std::endl;
        // std::cout << cursor << " " << pos << std::endl;
        // std::cout << "val: " << ((content[cursor] >> pos) & 1) << std::endl;
                
        return ((content[cursor] >> pos) & 1);
    }


    std::string to_string() const
    {
        std::stringstream result;
        result << std::setw(2*sizeof(BinWord)) << std::hex ;
        for(auto x_i : content)
            result << x_i << "  ";
        return result.str();
    }
    
    
    inline unsigned int get_msb() const
    {
        return msb;
    }
      

    inline bool is_zero() const
    {
        for (auto x_i : content)
            if (x_i != 0)
                return false;
        return true;
    }

    // !SUBSECTION! modifying the state

    
    inline void set_to_1(const unsigned int index)
    {
        if (index >= total_length)
            throw std::runtime_error("trying to set a bit with an index too high");
        else
        {
            unsigned int
                cursor = index / BLOCK_SIZE,
                pos    = index % BLOCK_SIZE;
            content[cursor] |= (1L << pos);
            if ((index > msb) || (msb == MSB_OF_ZERO))
                msb = index;
        }
    }

    inline void operator^=(const cpp_BigF2Vector & x)
    {
        if (total_length != x.size())
            throw std::runtime_error("xoring vectors of different size");
        else
        {
            for (unsigned int i=0; i<content.size(); i++)
                content[i] ^= x.content[i];
            set_msb();
        }
    }
};


// !SECTION! Operators and helper functions

inline cpp_BigF2Vector operator^(const cpp_BigF2Vector & x,
                                 const cpp_BigF2Vector & y)
{
    if (x.size() != y.size())
        throw std::runtime_error("xoring vectors of different size");
    else
    {
        std::vector<BinWord> result(x.content.begin(), x.content.end());
        for(unsigned int i=0; i<x.n_blocks(); i++)
            result[i] ^= y.content[i];
        return cpp_BigF2Vector(result, x.size());
    }
}


inline bool operator==(const cpp_BigF2Vector & x,
                       const cpp_BigF2Vector & y)
{
    return (x.content == y.content);
}


inline bool operator<(const cpp_BigF2Vector & x,
                      const cpp_BigF2Vector & y)
{
    return std::lexicographical_compare(
        x.content.begin(), x.content.end(),
        y.content.begin(), y.content.end()
        );
}

#endif
