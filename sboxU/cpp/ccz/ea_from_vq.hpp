#ifndef _CCZ_EA_FROM_VQ_
#define _CCZ_EA_FROM_VQ_


#include "../common.hpp"
#include "../core/include.hpp"
#include "./zeroes.hpp"
#include "./graph.hpp"
#include "../apn/include.hpp"
#include "../apn/ccz_class.hpp"

std::map<cpp_BinLinearBasis, std::vector<unsigned int>>cpp_image_of_space_by_group(std::vector<cpp_F2AffineMap> G, cpp_BinLinearBasis V);
std::vector<cpp_F2AffineMap> cpp_product_walsh_match(
    std::vector<cpp_F2AffineMap> G1,
    std::vector<cpp_F2AffineMap> G2,
    cpp_BinLinearBasis V1,
    cpp_BinLinearBasis V2);
bool cpp_product_walsh_match_any(
    const std::vector<cpp_F2AffineMap>& G1,
    const std::vector<cpp_F2AffineMap>& G2,
    const cpp_BinLinearBasis& V1,
    const cpp_BinLinearBasis& V2);

    
std::vector<cpp_F2AffineMap> cpp_ea_mapping_from_vq(
    const cpp_S_box f,
    const cpp_S_box g,
    const unsigned int n_threads,
    const std::string & mode = "standard");

#endif

