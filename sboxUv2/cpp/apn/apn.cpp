#include "apn.hpp"

cpp_S_box cpp_ortho_derivative(const cpp_S_box s)
{
    std::vector<BinWord> result(s.input_space_size(), 0);
    for (unsigned int a=1; a<s.input_space_size(); a++)
    {
        // getting the hyperplane
        std::vector<Integer> row(cpp_ddt_row(s, a));
        std::vector<Integer> hyperplane;
        hyperplane.reserve(s.size() / 2);
        for(unsigned int b=1; b<s.input_space_size(); b++)
            if (row[b] != row[0])
                hyperplane.push_back(b);

        // we return an empty list if the function is not APN, which
        // is equivalent to all rows having exactly half of their
        // elements be non-zero
        if (hyperplane.size() < (s.input_space_size()/2))
            return cpp_empty_S_box();

        // bruteforcing "ortho" until it is orthogonal to all elements
        // in the hyperplane
        BinWord ortho = 1;
        bool found = false;
        while ((not found) and (ortho < s.input_space_size()))
        {
            found = true;
            for(auto &b : hyperplane)
                if (cpp_scal_prod(ortho, b) == 0)
                {
                    found = false;
                    break;
                }
            if (not found)
                ortho += 1;
        }
        // if we couldn't find an element orthogonal to the
        // hyperplane, then it is not a hyperplane and we return an
        // empty function
        if (found)
            result[a] = ortho;
        else
            return cpp_empty_S_box();
    }
    return cpp_S_box(result);
}
