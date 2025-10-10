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
    graph.reserve(s.input_space_size());
    for(BinWord x=0; x<s.input_space_size(); x++)
        graph.push_back((x << n) | s[x]);
}


void cpp_FunctionGraph::set_content(
    const std::vector<BinWord> & v,
    BinWord _n
    )
{
    n = _n;
    mask = (1 << n) - 1;
    graph.assign(v.cbegin(), v.cend());
}


cpp_S_box cpp_FunctionGraph::get_S_box() const
{
    BinWord forbidden_value = (1 << n) + 1;
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
    BinWord forbidden_value = (1 << n) + 1;
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



// !SECTION! Helper functions 


cpp_BinLinearMap cpp_EA_mapping(
    const cpp_BinLinearMap &A,
    const cpp_BinLinearMap &B,
    const cpp_BinLinearMap &C
    )
{
    if (A.get_input_length() != C.get_input_length())
        throw std::runtime_error("in cpp_EA_BinLinearMap: incompatible block lengths");
    else
    {
        cpp_BinLinearMap
            A_inv = A.inverse(),
            C_prime = C * A_inv;
        std::vector<BinWord> images;
        for(unsigned int i=0; i<B.get_input_length(); i++)
        {
            images.push_back(B.image_vectors[i]);
        }
        for(unsigned int i=0; i<A_inv.get_input_length(); i++)
        {
            BinWord y = A_inv.image_vectors[i] << A_inv.get_output_length() ;
            y ^= C_prime.image_vectors[i];
            images.push_back(y);
        }
        return cpp_BinLinearMap(images);
    }
}

