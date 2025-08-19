#ifndef _CCZ_GRAPH_
#define _CCZ_GRAPH_

#include "../sboxU.hpp"
#include "../core/include.hpp"


class cpp_FunctionGraph
{
private:
    std::vector<BinWord> graph;
    BinWord mask;
    BinWord n;
    
public:
    cpp_FunctionGraph() : graph(), mask(0), n(0) {} ;

    cpp_FunctionGraph(const cpp_S_box & s);

    void init_from_S_box(const cpp_S_box & s);

    void set_content(const std::vector<BinWord> & v, BinWord _n);

    cpp_S_box get_S_box() const;

    cpp_S_box apply_basis_change(const std::vector<BinWord> &L) const;
};

#endif
