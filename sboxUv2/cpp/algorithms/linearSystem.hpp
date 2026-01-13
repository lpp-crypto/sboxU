#ifndef _LINEAR_SYSTEM_
#define _LINEAR_SYSTEM_

#include "../sboxU.hpp"
#include "../core/include.hpp"
#include "./bigvectors.hpp"


// !SECTION! Helper functions


bool maybe_add_vector(
    std::map<BinWord, cpp_BigF2Vector> & equations,
    cpp_BigF2Vector & new_eq
    );


bool surely_add_vector(
    std::vector<cpp_BigF2Vector> & equations,
    cpp_BigF2Vector & new_eq
    );


// !SECTION! The linear system itself 

class cpp_F2LinearSystem
{
private:
    BinWord n_var;
    bool echelonize;
    std::vector<cpp_BigF2Vector> all_equations;
    std::map<BinWord, cpp_BigF2Vector> echelonized_equations;
    std::map<BinWord, cpp_BigF2Vector> forbidden_solutions;
    
public:
    cpp_F2LinearSystem(
        const BinWord _n_var,
        const bool _echelonize
        ) :
        n_var(_n_var),
        echelonize(_echelonize)
    {}
    

    inline int rank() const
    {
        if (echelonize)
            return echelonized_equations.size();
        else
            return -1;
    }

    
    inline unsigned int size() const
    {
        if (echelonize)
            return echelonized_equations.size();
        else
            return all_equations.size();
    }
    
    
    inline bool add_equation(const std::vector<BinWord> & var_indices)
    {
        cpp_BigF2Vector new_eq(n_var);
        for (auto ind : var_indices)
        {
            new_eq.set_to_1(ind);
        }
        return add_equation(new_eq);
    }

    
    inline bool add_equation(cpp_BigF2Vector & eq)
    {
        if (echelonize)
            return maybe_add_vector(echelonized_equations, eq);
        else
        {
            all_equations.push_back(eq);
            return true;
        }
    }


    bool remove_solution(const std::vector<BinWord> & var_indices);

    bool remove_solution(cpp_BigF2Vector & eq)
    {
        return maybe_add_vector(forbidden_solutions, eq);
    }


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
    {
        ops.reserve(n_var * n_var);
    } 

    inline void push_back(
        const BinWord origin,
        const BinWord destination
        )
    {
        ops.emplace_back(std::make_pair(origin, destination));
    }

    
    inline cpp_BigF2Vector eval_canonical(const unsigned int i) const
    {
        cpp_BigF2Vector result(n_var);
        result.set_to_1(i);
        for (int k=ops.size()-1; k>=0; k--)
        {
            if (result.is_set(ops[k].second))
                result.content[BLOCK_INDEX(ops[k].first)] ^= ((BinWord)1 << BLOCK_POS(ops[k].first));
        }
        return result;
    }

    
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
