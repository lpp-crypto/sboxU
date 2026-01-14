#ifndef _LINEAR_SYSTEM_
#define _LINEAR_SYSTEM_

#include "../common.hpp"
#include "../core/include.hpp"
#include "./bigvectors.hpp"


// !SECTION! The linear system itself 

class cpp_F2LinearSystem
{
private:
    std::map<unsigned int, cpp_BigF2Vector> equations;
    std::map<unsigned int, cpp_BigF2Vector> forbidden_solutions;
    unsigned int n_var;
    
public:
    cpp_F2LinearSystem(const unsigned int _n_var) :
        n_var(_n_var) 
    {}
    
    inline unsigned int rank() const
    {
        return equations.size();
    }
    
    bool add_equation(const std::vector<unsigned int> & var_indices);

    bool remove_solution(const std::vector<unsigned int> & sol);
    
    std::vector<cpp_BigF2Vector> kernel();

    std::vector<Bytearray> kernel_as_bytes();

    std::string to_string() const;
    
    void init_image_vectors(
        std::vector<cpp_BigF2Vector> & vectors
        ) const;
};




// !SECTION! The XOR-sequence

class cpp_XorSequence
{
private:
    unsigned int n_var;
    std::vector<std::pair<unsigned int, unsigned int> > ops;
    
public:
    cpp_XorSequence(unsigned int _n_var) :
        n_var(_n_var),
        ops()
    {} ;

    inline void push_back(
        const unsigned int origin,
        const unsigned int destination
        )
    {
        ops.emplace_back(std::make_pair(origin, destination));
    };

    
    cpp_BigF2Vector eval_canonical(const unsigned int i) const
    {
        cpp_BigF2Vector result(n_var);
        result.set_to_1(i);
        for (int k=ops.size()-1; k>=0; k--)
        {
            if (result.is_set(ops[k].second))
            {
                unsigned int
                    cursor = ops[k].first / BLOCK_SIZE,
                    pos    = ops[k].first % BLOCK_SIZE;
                result.content[cursor] ^= ((BinWord)1 << pos);
            }
        }
        return result;
    };

    
    std::string to_string() const
    {
        std::stringstream result;
        result << "Xseq: " ;
        for (auto &p : ops)
            result << std::dec << "(" << p.first
                   << "," << p.second << ")  ";
        return result.str();
    }
};


#endif
