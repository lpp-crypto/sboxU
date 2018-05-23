/* Time-stamp: <2018-05-16 14:50:34 lperrin>
 *
 * LICENSE
 */ 



#include "sboxu_cpp.hpp"
using namespace boost::python;


// !SECTION! Wrapping functions to use python types 


// !SUBSECTION! Basic utils

bool is_permutation(list s)
{
    return is_permutation_cpp(lst_2_vec_BinWord(s));
}

list random_permutation(unsigned int n)
{
    return vec_2_lst_BinWord(random_permutation_cpp(n)); 
}


// !SUBSECTION! Linear Properties

list fourier_transform(const list& l)
{
    Sbox f(lst_2_vec_BinWord(l));
    std::vector<Integer> transform = walsh_spectrum_coord(f);
    return vec_2_lst_Integer(transform);
}

list invert_lat_fast(const list &t, const unsigned int n)
{
    std::vector<std::vector<Integer> > vec_t;
    for (unsigned int a=0; a<len(t); a++)
        vec_t.push_back(lst_2_vec_Integer(extract<list>(t[a])));
    return vec_2_lst_BinWord(invert_lat_cpp(vec_t, n)) ;
}


// !SUBSECTION! Linear/Affine equivalence

list linear_equivalence_fast(const list& l0, const list &l1)
{
    Sbox f(lst_2_vec_BinWord(l0)), g(lst_2_vec_BinWord(l1));
    std::vector<Sbox> mappings(linear_equivalence_cpp(f, g));
    list result;
    for (auto &m : mappings)
        result.append(vec_2_lst_BinWord(m));
    return result;
}

list le_class_representative(const list& l0)
{
    return vec_2_lst_BinWord(le_class_representative_cpp(lst_2_vec_BinWord(l0)));
}


// !SUBSECTION! CCZ

// list extract_vector_fast(const list& l, const Integer a) 
// {
//     return vec_2_lst_BinWord(extract_vector_cpp(lst_2_vec_BinWord(l), a)) ;
// }


list extract_bases_fast(const list& l,
                        const unsigned int dimension,
                        const unsigned int word_length,
                        const unsigned int n_threads) 
{
    std::vector<std::vector<BinWord> > bases = extract_bases_cpp(
        lst_2_vec_BinWord(l),
        dimension,
        word_length,
        n_threads);
    list result;
    for (auto &b : bases)
        result.append(vec_2_lst_BinWord(b));
    return result;
}


// !SECTION! Declaring all python-reachable function 

BOOST_PYTHON_MODULE(sboxu_cpp)
{
    // Utils
    def("oplus_cpp",
        oplus_cpp,
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
    def("parity",
        parity,
        args("x"),
        "Returns the parity of the input, i.e. its hamming weight modulo 2.");
    def("component",
        component,
        args("a", "f"),
        "Returns the Boolean function $x \\mapsto a.f(x)$, where the dot '.' denotes the scalar product in GF(2^n).");
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
    def("invert_lat_fast",
        invert_lat_fast,
        args("T", "n"),
        "Returns an S-Box s operating on n bits such that its LAT is the table T.");
    def("walsh_spectrum_fast",
        walsh_spectrum_fast,
        args("S", "n_threads"),
        "Returns a dictionnary d such that d[k] = #{(a, b), b != 0, k = \\sum_{x \\in \\{ 0,1 \\}^n} (-1)^{<a,x> + <b,S(x>)} which is computed using n_threads different threads.");
    def("fourier_transform",
        fourier_transform,
        args("f"),
        "Returns a list T such that T[a] = \\sum_{x \\in \\{ 0,1 \\}^n} (-1)^{<a,x> + S(x)}, where f is assumed to be a Boolean function.");
    def("lat_zeroes_fast",
        lat_zeroes_fast,
        args("S", "n", "n_threads"),
        "Returns all zeroes in the LAT of S as elements a||b of \ftwo^{2n}.");
    // Equivalence
    def("linear_equivalence_fast",
        linear_equivalence_fast,
        args("f", "g"),
        "Returns a list containing the LUT of two linear mappgings A and B such that f(x) = B(g(A(x))) or, if no such mappings exist, an empty list.");
    def("le_class_representative",
        le_class_representative,
        args("f"),
        "Returns the smallest representative of the linear-equivalence class of f.");
    def("extract_bases_fast",
        extract_bases_fast,
        args("V", "d", "n", "n_threads"),
        "Returns a list containing the minimal bases of all vector spaces of dimension d included in V.");
}

