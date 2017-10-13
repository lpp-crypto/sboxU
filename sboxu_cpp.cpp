/* Time-stamp: <2017-10-13 11:47:37 lperrin>
 *
 * LICENSE
 */ 



#include "sboxu_cpp.hpp"
using namespace boost::python;


// !SECTION! Wrapping functions to use python types 

// !SUBSECTION! Linear/Affine equivalence

list linear_equivalence_fast(const list& l0, const list &l1)
{
    Sbox f(lst_2_vec_int(l0)), g(lst_2_vec_int(l1));
    std::vector<Sbox> mappings(linear_equivalence_cpp(f, g));
    list result;
    for (auto &m : mappings)
        result.append(vec_2_lst_int(m)) ;
    return result;
}

list le_class_representative(const list& l0)
{
    return vec_2_lst_int(le_class_representative_cpp(lst_2_vec_int(l0)));
}


// !SECTION! Declaring all python-reachable function 

BOOST_PYTHON_MODULE(sboxu_cpp)
{
    // Utils
    def("oplus",
        oplus,
        args("x", "y"),
        "Returns the XOR of its (at most 32-bit) inputs.");
    def("hamming_weight",
        hamming_weight,
        args("x"),
        "Returns the Hamming weight of its input.");
    def("scal_prod",
        scal_prod,
        args("x", "y"),
        "Returns the scalar product in GF(2)^n of its inputs.");
    def("random_permutation",
        random_permutation,
        args("n"),
        "Returns a list containing the LUT of a permutation on n bits.");
    def("is_permutation",
        is_permutation,
        args("S"),
        "Returns True if and only if the list S corresponds to a permutation.");
    // Differential properties
    def("ddt",
        ddt,
        args("S"),
        "Returns a list of lists d such that d[a][b] = #{x, S(x ^ a) ^ S(x) = b}");
    def("differential_spectrum_fast",
        differential_spectrum_fast,
        args("S", "n_threads"),
        "Returns a dictionnary d such that d[k] = #{(a, b), a != 0, S(x ^ a) ^ S(x) = b has k solutions} which is computed using n_threads different threads.");
    // Linear properties
    def("lat",
        lat,
        args("S"),
        "Returns a list of lists l such that l[a][b] = \\sum_{x \\in \\{ 0,1 \\}^n} (-1)^{<a,x> + <b,S(x>)}.");
    def("walsh_spectrum_fast",
        walsh_spectrum_fast,
        args("S", "n_threads"),
        "Returns a dictionnary d such that d[k] = #{(a, b), b != 0, k = \\sum_{x \\in \\{ 0,1 \\}^n} (-1)^{<a,x> + <b,S(x>)} which is computed using n_threads different threads.");
    // Equivalence
    def("linear_equivalence_fast",
        linear_equivalence_fast,
        args("f", "g"),
        "Returns a list containing the LUT of two linear mappgings A and B such that f(x) = B(g(A(x))) or, if no such mappings exist, an empty list.");
    def("le_class_representative",
        le_class_representative,
        args("f"),
        "Returns the smallest representative of the linear-equivalence class of f.");
}

