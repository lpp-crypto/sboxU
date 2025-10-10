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

    cpp_S_box get_ccz_equivalent_function(const cpp_BinLinearMap &L) const;
};


/* Returns a BinLinearMap corresponding to the n+m bit linear application that must be applied to the graph of a function F in order to obtain the graph of the function BoFoA + C.
 */
cpp_BinLinearMap cpp_EA_mapping(
    const cpp_BinLinearMap &A,
    const cpp_BinLinearMap &B,
    const cpp_BinLinearMap &C
    );


#endif
