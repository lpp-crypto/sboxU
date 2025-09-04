#ifndef _CCZ_LINEAR_EQUIV_
#define _CCZ_LINEAR_EQUIV_

#include "../sboxU.hpp"
#include "../core/include.hpp"



/** A and B are reset to contain the linear mappings such that repr = BoFoA.
 *
 * A natural approach would have involved the output of a std::tuple, but these cannot be interfaced with cython.
 */
cpp_S_box cpp_le_class_representative(
    const cpp_S_box f,
    cpp_BinLinearMap & A,
    cpp_BinLinearMap & B
    );


cpp_S_box cpp_le_class_representative(
    const cpp_S_box f
    )

{
    cpp_BinLinearMap A, B;
    return cpp_le_class_representative(f, A, B);
};


#endif /* _SBOXU_CPP_EQUIV_H_ */
