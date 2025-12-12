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

    // !TODO! add logic to simplify this list using the same tricks as for equation addition 
    if (sol.size() > 0)
    {
        cpp_BigF2Vector f(n_var);
        for (auto ind : sol)
        {
            f.set_to_1(ind);
        }
        forbidden_solutions.push_back(f);
    }
    // if (sol.size() == 0)
    //     return;
    // else
    // {
    //     // if ((sol.size() % 2) == 1)
    //     // {
    //     //     // in this case, the scalar product sol*sol is 1, meaning that
    //     //     // adding sol as a row of the system will prevent it from
    //     //     // being a solution
    //     //     add_equation(sol);
    //     // }
    //     else
    //     {
    //         // otherwise, we remove one entry
    //         add_equation(std::vector<unsigned int>(sol.begin(), sol.end()-1));
    //     }
    // }
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
        for (unsigned int i=0; i<vectors.size(); i++)
            if (vectors[i].is_non_zero())
                for (unsigned int j=i+1; j<vectors.size(); j++)
                {
                    cpp_BigF2Vector x = vectors[i] ^ vectors[j];
                    if (x < vectors[j])
                    {
                        vectors[j] = x;
                        seq.push_back(i,j);
                    }
                }
        // and then we rebuild the kernel vectors
        std::vector<cpp_BigF2Vector> ker;
        for (unsigned int z=0; z<vectors.size(); z++)
            if (vectors[z].is_zero())            
                ker.push_back(seq.eval_canonical(z));
        // finally, we remove the contribution of the forbidden solutions
        unsigned int n_added = forbidden_solutions.size();
        if (n_added > 0)
        {
            std::vector<cpp_BigF2Vector> ker_prime(forbidden_solutions.begin(),
                                                   forbidden_solutions.end());
            ker_prime.insert(ker_prime.end(), ker.begin(), ker.end());
            ker.clear();
            for (unsigned int i=0; i<ker_prime.size(); i++)
                if (ker_prime[i].is_non_zero())
                {
                    if (i >= n_added)
                        ker.push_back(ker_prime[i]);
                    for (unsigned int j=i+1; j<ker_prime.size(); j++)
                    {
                        cpp_BigF2Vector x = ker_prime[i] ^ ker_prime[j];
                        if (x < ker_prime[j])
                            ker_prime[j] = x;
                    }
                }
        }
        // we finally build the result
        return ker;
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
    unsigned int cursor = 0, pos=0;
    BinWord mask = 1;
    for (auto eq : equations)
    {
        for (unsigned int i=0; i<n_var; i++)
            if (eq.second.is_set(i))
                vectors[i].content[cursor] |= mask;
        pos ++;
        mask <<= 1;
        if (pos == BLOCK_SIZE)
        {
            mask = 1;
            cursor ++;
            pos = 0;
        }
    }
    for(unsigned int i=0; i<vectors.size(); i++)
        vectors[i].set_msb();
}
