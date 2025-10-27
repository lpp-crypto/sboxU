#include "./graph.hpp"


// !SECTION! The cpp_FunctionGraph class 

cpp_FunctionGraph::cpp_FunctionGraph(const cpp_S_box & s) :
    graph(),
    mask(0),
    n(0)
{
    init_from_S_box(s);
}


void cpp_FunctionGraph::init_from_S_box(
    const cpp_S_box & s
    )
{
    n = s.get_output_length();
    mask = (1 << n) - 1;
    graph.assign(s.input_space_size(), 0);
    for(BinWord x=0; x<s.input_space_size(); x++)
        graph[x] = (x << n) | s[x];
    // no need to sort the graph: it is sorted by construction
}


void cpp_FunctionGraph::set_content(
    const std::vector<BinWord> & v,
    BinWord _n
    )
{
    n = _n;
    mask = (1 << n) - 1;
    graph.assign(v.cbegin(), v.cend());
    // we need to make sure that the graph is sorted
    std::sort(graph.begin(), graph.end());
}


cpp_S_box cpp_FunctionGraph::get_cpp_S_box() const
{
    BinWord forbidden_value = 1 << (n + 1);
    std::vector<BinWord> lut(graph.size(), forbidden_value);
    for (auto & entry : graph)
    {
        BinWord
            x = entry >> n,
            y = entry & mask;
        if (lut[x] != forbidden_value)
            // problem: this entry has already been set
            return cpp_empty_S_box();
        else
            lut[x] = y;
    }
    return cpp_S_box(lut);
}


cpp_S_box cpp_FunctionGraph::get_ccz_equivalent_function(
    const cpp_BinLinearMap &L
    ) const
{
    BinWord forbidden_value = 1 << (n + 1);
    std::vector<BinWord> lut(graph.size(), forbidden_value);
    for (auto & entry : graph)
    {
        BinWord
            z = L(entry),
            x = z >> n,
            y = z & mask;
        if (lut[x] != forbidden_value)
            // problem: this entry has already been set
            return cpp_empty_S_box();
        else
            lut[x] = y;
    }
    return cpp_S_box(lut);
}


std::vector<BinWord> cpp_FunctionGraph::xor_equivalence(const cpp_FunctionGraph & G) const
{
    std::vector<BinWord> result;
    if (graph.size() != G.size())
        return result;
    for (auto p_0 : graph)
    {
        bool valid_offset = true;
        BinWord offset = p_0 ^ G[0];
        for (auto io_pair : graph)
            if (not G.contains(io_pair ^ offset))
            {
                valid_offset = false;
                break;
            }
        if (valid_offset)
            result.push_back( offset );
    }       
    return result;
}


cpp_FunctionGraph cpp_FunctionGraph::image_by(const cpp_BinLinearMap & L) const
{
    BinWord forbidden_value = 1 << (n + 1);
    std::vector<BinWord> new_graph;
    new_graph.reserve(graph.size());
    for(auto & entry : graph)
        new_graph.push_back(L(entry));
    cpp_FunctionGraph result;
    result.set_content(new_graph, n);
    return result;
}
    

// !SECTION! Helper functions 


cpp_BinLinearMap cpp_EA_mapping(
    const cpp_BinLinearMap &A,
    const cpp_BinLinearMap &B,
    const cpp_BinLinearMap &C
    )
{
    if (A.get_input_length() != A.get_output_length()) 
        throw std::runtime_error("in cpp_EA_BinLinearMap: A cannot be invertible");
    else if (B.get_input_length() != B.get_output_length()) 
        throw std::runtime_error("in cpp_EA_BinLinearMap: B cannot be invertible");
    else if (A.get_input_length() < C.get_input_length())
        throw std::runtime_error("in cpp_EA_BinLinearMap: A and C have incompatible input length");
    else if (B.get_output_length() < C.get_output_length())
        throw std::runtime_error("in cpp_EA_BinLinearMap: B and C have incompatible output length");
    else
    {
        cpp_BinLinearMap
            A_inv = A.inverse(),
            C_prime = C * A_inv;
        std::vector<BinWord> images;
        for(unsigned int i=0; i<B.get_input_length(); i++)
            images.push_back(B(1 << i));
        for(unsigned int i=0; i<A.get_output_length(); i++)
        {
            BinWord x = 1 << i;
            images.push_back((A_inv(x) << B.get_output_length()) | C_prime(x));
        }
        return cpp_BinLinearMap(images);
    }
}

// Returns the block decomposition of a 2n x 2n matrix into four n x n blocks.
std::vector<cpp_BinLinearMap> ccz_block_decomposition(const cpp_BinLinearMap &L) {
        if(L.get_input_length() != L.get_output_length() || L.get_input_length() & ((Integer) 1))
            return {};

        // A D
        // C B
        std::vector<std::vector<BinWord>> images_of_each_block(4); // A B C D
        Integer n = L.get_input_length() >> 1;

        BinWord v, w;
        for (int i = 0; i < n; i++) {
            v = L.image_vectors[i];
            w = L.image_vectors[(n + i)];
            images_of_each_block[0].push_back((v << n) >> n); //A
            images_of_each_block[2].push_back(v >> n); // C

            images_of_each_block[3].push_back((w << n) >> n); // D
            images_of_each_block[1].push_back(w >> n); // B
        }

        std::vector<cpp_BinLinearMap> block_decomposition;
        for (int i = 0; i < images_of_each_block.size(); i++)
            block_decomposition.push_back(cpp_BinLinearMap(images_of_each_block[i], n, n));
        return block_decomposition;
}


