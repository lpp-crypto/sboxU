#ifndef _LINEAR_SYSTEM_
#define _LINEAR_SYSTEM_

#include "../sboxU.hpp"
#include "../core/include.hpp"
#include "./bigvectors.hpp"


// !SECTION! The linear system itself 

class cpp_F2LinearSystem
{
private:
    std::map<BinWord, cpp_BigF2Vector> equations;
    std::map<BinWord, cpp_BigF2Vector> forbidden_solutions;
    BinWord n_var;
    
public:
    cpp_F2LinearSystem(const BinWord _n_var) :
        n_var(_n_var) 
    {}
    
    inline BinWord rank() const
    {
        return equations.size();
    }
    
    bool add_equation(const std::vector<BinWord> & var_indices);
    
    bool add_equation(cpp_BigF2Vector & eq);

    bool remove_solution(const std::vector<BinWord> & var_indices);

    bool remove_solution(cpp_BigF2Vector & eq);

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
    BinWord n_var;
    std::vector<std::pair<BinWord, BinWord> > ops;
    
public:
    cpp_XorSequence(BinWord _n_var) :
        n_var(_n_var),
        ops()
    {} ;

    inline void push_back(
        const BinWord origin,
        const BinWord destination
        )
    {
        ops.emplace_back(std::make_pair(origin, destination));
    };

    
    inline cpp_BigF2Vector eval_canonical(const unsigned int i) const
    {
        cpp_BigF2Vector result(n_var);
        result.set_to_1(i);
        for (int k=ops.size()-1; k>=0; k--)
        {
            if (result.is_set(ops[k].second))
            {
                BinWord
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
