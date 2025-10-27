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

    cpp_S_box get_cpp_S_box() const;

    cpp_S_box get_ccz_equivalent_function(const cpp_BinLinearMap &L) const;

    std::vector<BinWord> xor_equivalence(const cpp_FunctionGraph & G) const;

    cpp_FunctionGraph image_by(const cpp_BinLinearMap & L) const;

    inline bool contains(const BinWord x) const
    {
        return std::binary_search(graph.begin(), graph.end(), x);
    }

    inline bool operator==(const cpp_FunctionGraph & g) const
    {
        BinWord x = 0;
        for(auto & entry : g)
        {
            if (graph[x] != entry)
                return false;
            x ++;
        }
        return true;
    }       

    inline BinWord operator[](const BinWord & x) const
    {
        return graph[x];
    }
    
    inline Integer size() const
    {
        return graph.size();
    }
    
    inline std::vector<BinWord>::const_iterator begin() const
    {
        return graph.cbegin();
    }
    
    inline std::vector<BinWord>::const_iterator end() const
    {
        return graph.cend();
    }

};


/* Returns a BinLinearMap corresponding to the n+m bit linear application that must be applied to the graph of a function F in order to obtain the graph of the function BoFoA + C.
 */
cpp_BinLinearMap cpp_EA_mapping(
    const cpp_BinLinearMap &A,
    const cpp_BinLinearMap &B,
    const cpp_BinLinearMap &C
    );


#endif
