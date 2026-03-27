#include "./ortho_derivative.hpp"

#include <random>
#include <unordered_set>

// !SECTION! The ortho-derivative itself

cpp_S_box cpp_ortho_derivative(const cpp_S_box &s)
{
    const Integer n = s.get_input_length();
    const Integer N = s.input_space_size();

    std::vector<BinWord> result(N, 0);

    std::mt19937_64 rng(std::random_device{}());
    std::uniform_int_distribution<BinWord> dist(0, N - 1);

    for (BinWord a = 1; a < N; a++)
    {
        std::unordered_set<BinWord> image_values;
        cpp_BinLinearBasis basis;
        BinWord y0 = 0; // Some fixed point that is contained in the underlying vec space of Im(D_aF)
        bool found_y0 = false;

        while ((basis.rank() < n - 1))
        {
            BinWord x = dist(rng);
            BinWord y = s[x ^ a] ^ s[x];
            if (!found_y0)
            {
                y0 = y;
                found_y0 = true;
                image_values.insert(y);
                continue;
            }
            if (!image_values.insert(y).second) // that is, if y is already in image_values
                continue;
            basis.add_to_span(y ^ y0); // Note that we already know y \neq y_0.
        }
        // Build the linear system <u, v> = 0 for all basis vectors v.
        cpp_F2LinearSystem eqs(n, true);
        for (const BinWord v : basis.get_basis())
        {
            std::vector<BinWord> positions;
            for (Integer i = 0; i < n; i++)
                if ((v >> i) & 1)
                    positions.push_back(i);
            eqs.add_equation(positions);
        }

        std::vector<cpp_BigF2Vector> ker = eqs.kernel();
        if (ker.size() != 1)
            return cpp_empty_S_box(); // Something went wrong or given function was not crooked.

        BinWord u = 0;
        for (Integer i = 0; i < n; i++)
            if (ker[0].is_set(i))
                u |= ((BinWord)1 << i);
        if (u == 0)
            return cpp_empty_S_box(); // Should only occur when a = 0.

        result[a] = u;
    }

    return cpp_S_box(result);
}

// !SECTION! Ortho-integration

cpp_S_box cpp_S_box_from_bits(const cpp_BigF2Vector &v)
{
    unsigned int n = 1;
    while ((n * (1 << n) != v.size()) && (n < 32))
        n++;
    if (n == 32)
        throw std::runtime_error("input vector for cpp_S_box_from_bits has the wrong length");
    else
    {
        std::vector<BinWord> result(1 << n, 0);
        for (BinWord x = 0; x < result.size(); x++)
            for (unsigned int i = 0; i < n; i++)
            {
                result[x] <<= 1;
                if (v.is_set(i + x * n))
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

cpp_S_box cpp_ortho_integral(const cpp_S_box &s)
{
    BinWord n = s.get_input_length();
    // the system is expected to be big so we echelonize it to save space
    cpp_F2LinearSystem eqs(n * (1 << n), true);
    // ensuring that the solution F is such that F(0)=0
    for (BinWord i = 0; i < n; i++)
        eqs.add_equation(std::vector<BinWord>(1, i));
    // removing linear solutions
    for (BinWord i = 0; i < n; i++)
        for (BinWord j = 0; j < n; j++)
        {
            cpp_BigF2Vector L_ij(n * (1 << n));
            for (BinWord x = 0; x < s.input_space_size(); x++)
                if (((x >> i) & 1) == 1)
                    L_ij.set_to_1(j + x * n);
            eqs.remove_solution(L_ij);
        }
    // adding main equations
    for (BinWord a = 1; a < s.input_space_size(); a++)
    {
        for (BinWord x = 0; x < s.input_space_size(); x++)
            if (x != a)
            {
                std::vector<BinWord> positions;
                for (BinWord i = 0; i < n; i++)
                    if (((s[x] >> i) & 1) == 1)
                    {
                        positions.push_back(i + n * cpp_oplus(x, a)); // F_i(x+a)
                        positions.push_back(i + n * x);               // F_i(x)
                        positions.push_back(i + n * a);               // F_i(a)
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
