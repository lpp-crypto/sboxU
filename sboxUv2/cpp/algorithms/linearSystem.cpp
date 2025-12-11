#include "./linearSystem.hpp"


// !SECTION! Implementation of the cpp_F2LinearSystem class

// !SUBSECTION! Ensuring that only useful equations are kept

bool cpp_F2LinearSystem::add_equation(
    const std::vector<unsigned int> & var_indices
    )
{
    // building new equation
    cpp_BigF2Vector new_eq(n_var);
    for (auto ind : var_indices)
    {
        new_eq.set_to_1(ind);
    }
    // adding the equation as a horizontal vector: we perform a
    // reduction on the rows in the process, which will *not* change the
    // structure of the right kernel of the matrix.
    // -- reduction step
    for (auto  eq : equations)
    {        
        if (new_eq.get_msb() > eq.first)
        {
            // we reduce the current equation by the previous one
            if (new_eq.is_set(eq.first))
                new_eq ^= eq.second;
        }
        else if (new_eq.get_msb() == eq.first)
        {
            if (new_eq == eq.second)
            {
                // the new equation was already in the span: we stop
                // trying to add it
                return false;
            }
            else
            {
                cpp_BigF2Vector x = new_eq ^ eq.second;
                if (x < new_eq)
                    new_eq = x;
            }
        }
        else
            break;
    }
    // -- adding the equation
    unsigned int m = new_eq.get_msb();
    equations[m] = new_eq;
    // -- reducing the remaining equations
    for (auto eq : equations)
        if (eq.first > new_eq.get_msb())
        {
            // we reduce the old equation by the new one if the MSB of
            // the new equation is active in it
            if (eq.second.is_set(m))
                equations[eq.first] ^= new_eq;
        }
    return true;
}


// !SUBSECTION! Removing unwanted solutions

void cpp_F2LinearSystem::remove_solution(
    const std::vector<unsigned int> & sol
    )
{
    if (sol.size() == 0)
        return;
    else
    {
        if ((sol.size() % 2) == 1)
        {
            // in this case, the scalar product sol*sol is 1, meaning that
            // adding sol as a row of the system will prevent it from
            // being a solution
            add_equation(sol);
        }
        else
        {
        // otherwise, we remove one entry
            add_equation(std::vector<unsigned int>(sol.begin()+1, sol.end()));
        }
    }
}


// !SUBSECTION! Computing the kernel 

std::vector<cpp_BigF2Vector> cpp_F2LinearSystem::kernel()
{
    std::vector<cpp_BigF2Vector> result;
    if (equations.size() == n_var)
        // the equations are linearly independent by construction, so
        // if there are as many as variables then the kernel is
        // trivial
        return result;
    else
    {
        // we start by preparing the system to solve by computing the
        // column vectors
        std::vector<cpp_BigF2Vector> vectors;
        init_image_vectors(vectors);
        // we then move on to the elimination itself
        cpp_XorSequence seq(n_var);
        std::vector<unsigned int> zero_preimages;
        zero_preimages.reserve(n_var - equations.size());
        for (unsigned int i=0; i<vectors.size(); i++)
        {
            if (vectors[i].is_zero())
                zero_preimages.push_back(i);
            else
                for (unsigned int j=i+1; j<vectors.size(); j++)
                {
                    if ((vectors[i].get_msb() <= vectors[j].get_msb())
                        &&
                        (vectors[j].is_set(vectors[i].get_msb())))
                    {
                        // case where we reduce vectors[j]
                        vectors[j] ^= vectors[i];
                        seq.push_back(i, j);
                    }
                    else if ((vectors[i].get_msb() > vectors[j].get_msb())
                             &&
                             (vectors[i].is_set(vectors[j].get_msb())))
                    {
                        // case where we reduce vectors[i]
                        vectors[i] ^= vectors[j];
                        seq.push_back(j, i);
                    }
                }
        }   
        // and then we rebuild the kernel vectors
        result.reserve(n_var - equations.size());
        for (auto & z : zero_preimages)
        {
            result.push_back(seq.eval_canonical(z));
        }
        return result;
    }
} 


std::vector<Bytearray> cpp_F2LinearSystem::kernel_as_bytes()
{
    std::vector<cpp_BigF2Vector> ker = kernel();
    std::vector<Bytearray> result;
    for (auto v : ker)
        result.push_back(v.to_bytes());
    return result;
}


// !SUBSECTION! Helpers 

std::string cpp_F2LinearSystem::to_string() const
{
    std::stringstream result;
    for(auto & eq : equations)
    {
        result << std::setw(5) << std::dec << eq.first << " | ";
        for (unsigned int i=0; i<n_var; i++)
            if (eq.second.is_set(i))
                result << "1 ";
            else
                result << "0 ";
        result << std::endl;
    }
    return result.str();
}


void cpp_F2LinearSystem::init_image_vectors(
    std::vector<cpp_BigF2Vector> & vectors
    ) const
{
    vectors = std::vector<cpp_BigF2Vector> (
        n_var,
        cpp_BigF2Vector(equations.size())
        );
    unsigned int cursor = 0;
    BinWord mask = 1;
    for (auto eq : equations)
    {
        for (unsigned int i=0; i<n_var; i++)
            if (eq.second.is_set(i))
                vectors[i].content[cursor] |= mask;
        mask <<= 1;
        if (mask == 0)
        {
            mask = 1;
            cursor ++;
        }
    }
    for(unsigned int i=0; i<vectors.size(); i++)
        vectors[i].set_msb();
}
