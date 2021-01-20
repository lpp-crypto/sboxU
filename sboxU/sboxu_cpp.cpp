/* Time-stamp: <2021-01-19 17:01:53 leo>
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

dict cycle_type(list s)
{
    std::vector<Integer> lut(lst_2_vec_Integer(s));
    Integer i = 0;
    std::map<Integer, Integer> result;
    do
    {
        if (lut[i] != -1)
        {
            unsigned int length = 0;
            Integer j = i, prev_j = 0;
            do
            {
                prev_j = j;
                j = lut[j];
                lut[prev_j] = -1;
                length += 1;
            } while (j != i);
            result[length] += 1;
        }
        i ++;
    } while (i < lut.size());    
    dict result_py;
    for (auto it : result)
        result_py[it.first] = it.second;
    return result_py;
}


// !SUBSECTION! Linear Properties

list fourier_transform(const list& l)
{
    Sbox f(lst_2_vec_BinWord(l));
    std::vector<Integer> transform = walsh_spectrum_coord(f);
    return vec_2_lst_Integer(transform);
}

list invert_lat_fast(const list &t,
                     const unsigned int n)
{
    std::vector<std::vector<Integer> > vec_t;
    for (unsigned int a=0; a<len(t); a++)
        vec_t.push_back(lst_2_vec_Integer(extract<list>(t[a])));
    return vec_2_lst_BinWord(invert_lat_cpp(vec_t, n)) ;
}

list lat_zeroes_fast(const list &l,
                     const unsigned int n,
                     const unsigned int n_threads)

{
    std::vector<BinWord> s(lst_2_vec_BinWord(l)) ;
    return vec_2_lst_BinWord(lat_zeroes_cpp(s, n, n_threads));
}

list projected_lat_zeroes_fast(const list &l,
                               const unsigned int n_threads)

{
    std::vector<BinWord> s(lst_2_vec_BinWord(l)) ;
    return vec_2_lst_BinWord(projected_lat_zeroes_cpp(s, n_threads));
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

list extract_vector(const list& l, const BinWord a) 
{
    return vec_2_lst_BinWord(extract_vector_cpp(lst_2_vec_BinWord(l), a)) ;
}


list extract_bases_fast(const list& l,
                        const unsigned int dimension,
                        const unsigned int word_length,
                        unsigned int n_threads,
                        const str end_condition)
{
    std::vector<BinWord> space(lst_2_vec_BinWord(l)) ;
    std::vector<std::vector<BinWord> > bases = extract_bases_cpp(
        space,
        dimension,
        word_length,
        n_threads,
        extract<std::string>(end_condition));
    list result;
    for (auto &b : bases)
        result.append(vec_2_lst_BinWord(b));
    return result;
}


list extract_affine_bases_fast(const list& l,
                               const unsigned int dimension,
                               const unsigned int word_length,
                               unsigned int n_threads,
                               const str end_condition) 
{
    std::vector<BinWord> space(lst_2_vec_BinWord(l)) ;
    std::vector<std::vector<BinWord> > bases = extract_affine_bases_cpp(
        space,
        dimension,
        word_length,
        n_threads,
        extract<std::string>(end_condition));
    list result;
    for (auto &b : bases)
        result.append(vec_2_lst_BinWord(b));
    return result;
}


list get_lat_zeroes_spaces_fast(const list& l,
                                const unsigned int n,
                                const unsigned int n_threads)
{
    
    std::vector<BinWord> s(lst_2_vec_BinWord(l));
    std::vector<BinWord> zeroes(lat_zeroes_cpp(s, n, n_threads)) ;
    std::vector<std::vector<BinWord> > bases = extract_bases_cpp(
        zeroes,
        n,
        2*n,
        n_threads,
        std::string("fixed dimension"));
    list result;
    for (auto &b : bases)
        result.append(vec_2_lst_BinWord(b));
    return result;
}


// !SECTION! Other simple functions

Integer rank_of_vector_set_cpp(const list& V)
{
    Integer result = 0;
    std::vector<BinWord> l(lst_2_vec_BinWord(V));
    for (unsigned int i=0; i<l.size(); i++)
    {
        if (l[i] > 0)
        {
            result ++;
            for (unsigned int j=i+1; j<l.size(); j++)
            {
                BinWord y = l[i] ^ l[j];
                if (y < l[j])
                    l[j] = y;
            }
        }
    }
    return result;
}


bool rank_deficit_of_vector_set_is_at_most_cpp(const list& V, const Integer target)
{
    Integer count = 0;
    std::vector<BinWord> l(lst_2_vec_BinWord(V));
    for (unsigned int i=0; i<l.size(); i++)
    {
        if (l[i] == 0)
        {
            count ++;
            if (count > target)
                return false;
        }
        else
        {
            for (unsigned int j=i+1; j<l.size(); j++)
            {
                BinWord y = l[i] ^ l[j];
                if (y < l[j])
                    l[j] = y;
            }
        }
    }
    return true;
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
    def("cycle_type",
        cycle_type,
        args("S"),
        "If `S` is the LUT of a permutation P, returns a dictionnary `d` such that `d[k]` = l if and only if the permutation P has exactly `l` cycles of length `k`.");

// Differential properties
    
    def("ddt",
        ddt,
        args("S"),
        "Returns a list of lists d such that d[a][b] = #{x, S(x ^ a) ^ S(x) = b}");
    def("differential_spectrum_fast",
        differential_spectrum_fast,
        args("S", "n_threads"),
        "Returns a dictionnary d such that d[k] = #{(a, b), a != 0, S(x ^ a) ^ S(x) = b has k solutions} which is computed using n_threads different threads.");

    def("ortho_derivative",
        ortho_derivative,
        args("S"),
        "Returns a list containing the LUT of the ortho-derivative of S if S is both crooked and APN, an empty list otherwise.");

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
    def("projected_lat_zeroes_fast",
        projected_lat_zeroes_fast,
        args("S", "n_threads"),
        "Returns the projected of the set of all zeroes in the LAT of S, i.e. all a such that LAT[a,b]==0 for some b.");

// Boomerang properties
    
    def("bct",
        bct,
        args("S"),
        "Returns a list of lists b such that b[a][b] = #{x, S^-1(S(x)^b) ^ S^-1(S(x^a)^b) = a}");
    def("bct_spectrum_fast",
        bct_spectrum_fast,
        args("S", "n_threads"),
        "Returns a dictionnary b such that b[k] = #{(a, b), a != 0, b != 0, S^-1(S(x)^b) ^ S^-1(S(x^a)^b) = a} has k solutions} which is computed using n_threads different threads.");

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
    def("extract_affine_bases_fast",
        extract_affine_bases_fast,
        args("V", "d", "n", "n_threads"),
        "Returns a list containing the canonical bases of all affine spaces of dimension d included in V.");
    def("get_lat_zeroes_spaces_fast",
        get_lat_zeroes_spaces_fast,
        args("S", "n", "n_threads"),
        "Returns a list containing the minimal bases of all vector spaces of dimension n included in the Walsh zeroes of S.");

// Other functions
    def("extract_vector",
        extract_vector,
        args("l", "a"),
        "Extracts the vector a from the list l in the sense of [BonPerTia19]."),
    def("rank_of_vector_set_cpp",
        rank_of_vector_set_cpp,
        args("V"),
        "Returns the rank of the binary representations of the unsigned integers in V.");
    def("rank_deficit_of_vector_set_is_at_most_cpp",
        rank_deficit_of_vector_set_is_at_most_cpp,
        args("V", "target"),
        "Returns true if and only if the c-r>=target, where `r` is rank of the binary representations of the unsigned integers in V and `c` is the number of elements in `V`.");
}

