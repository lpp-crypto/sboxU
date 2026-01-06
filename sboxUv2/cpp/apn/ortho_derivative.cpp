#include "./ortho_derivative.hpp"


// !SECTION! The ortho-derivative itself

cpp_S_box cpp_ortho_derivative(const cpp_S_box &s)
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



// !SECTION! Ortho-integration

cpp_S_box cpp_S_box_from_bits(const cpp_BigF2Vector & v)
{
    unsigned int n = 1;
    while ((n*(1 << n) != v.size()) && (n < 32))
        n ++;
    if (n == 32)
        throw std::runtime_error("input vector for cpp_S_box_from_bits has the wrong length");
    else
    {
        std::vector<BinWord> result(1 << n, 0);
        for (BinWord x=0; x<result.size(); x++)
            for (unsigned int i=0; i<n; i++)
            {
                result[x] <<= 1;
                if (v.is_set(i + x*n))
                    result[x] |= 1;
            }
        return cpp_S_box(result, n, n);
    }
}


// cpp_BigF2Vector cpp_bits_from_S_box(const cpp_S_box & s)
// {
//     unsigned int n = s.get_input_length();
//     if (s.get_output_length() != n)
//         throw std::runtime_error("ortho_integral cannot be computed when input and output length do not match.");
//     cpp_BigF2Vector result(n * (1 << n));
//     for (BinWord x=0; x<result.size(); x++)
//         for (unsigned int i=0; i<n; i++)
//             if (((s[x] >> i) & 1) == 1)
//                 result.set_to_1(i);
//     return result;
// }


cpp_S_box cpp_ortho_integral(const cpp_S_box & s)
{
    unsigned int n = s.get_input_length();
    cpp_F2LinearSystem eqs(n*(1 << n));
    // ensuring that the solution F is such that F(0)=0
    for(unsigned int i=0; i<n; i++)
        eqs.add_equation(std::vector<unsigned int>(1, i));
    // removing linear solutions
    for(unsigned int i=0; i<n; i++)
        for(unsigned int j=0; j<n; j++)
        {
            cpp_BigF2Vector L_ij(n*(1 << n));
            for (BinWord x=0; x<s.input_space_size(); x++)
                if (((x >> i) & 1) == 1)
                    L_ij.set_to_1(j + x*n);
            eqs.remove_solution(L_ij);
        }
    // adding main equations
    for(BinWord a=1; a<s.input_space_size(); a++)
    {
        for(BinWord x=0; x<s.input_space_size(); x++)
            if (x != a)
            {
                std::vector<unsigned int> positions;
                for(unsigned int i=0; i<n; i++)
                    if (((s[x] >> i) & 1) == 1)
                    {
                        positions.push_back(i + n*cpp_oplus(x, a)); // F_i(x+a)
                        positions.push_back(i + n*x);               // F_i(x)
                        positions.push_back(i + n*a);               // F_i(a)
                    }
                eqs.add_equation(positions);
            }
    }
    // solving and outputting
    std::vector<cpp_BigF2Vector> ker = eqs.kernel();
    if (ker.size() > 1)
    {
        std::cout << "[BIG RESULT]\n" 
                  << " the function with the following LUT admits multiple ("
                  << ker.size()
                  << ") non-trivial ortho-integrals!\n"
                  << s.content_string_repr()
                  << std::endl;
        return cpp_S_box_from_bits(ker[0]);
    }
    else if (ker.size() == 0)
        return cpp_empty_S_box();
    else
        return cpp_S_box_from_bits(ker[0]);
}
