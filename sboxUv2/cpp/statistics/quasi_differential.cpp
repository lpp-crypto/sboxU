#include "quasi_differential.hpp"

Integer cpp_qddt_coeff(const cpp_S_box &s, const BinWord a, const BinWord b, const BinWord u, const BinWord v){
    std::vector<BinWord> xddt_a_b=cpp_xddt_entry(s,a,b); 
    Integer res = 0;
    for (const BinWord &x : xddt_a_b ){
        if (cpp_scal_prod(u,x)^cpp_scal_prod(v,s[x])==0){
            res ++;
        }
        else {
            res --;
        }
    }
    return res;
}

std::vector<std::vector<Integer>> cpp_qddt_fixed_differential(const cpp_S_box &s, const BinWord a, const BinWord b){
    std::vector<BinWord> xddt_a_b=cpp_xddt_entry(s,a,b); 
    std::vector< std::vector<Integer>> result;
    for (BinWord u = 0; u< s.input_space_size(); u++){
        std::vector<Integer> row;
        for (BinWord v = 0; v< s.output_space_size(); v++){ // Duplicating code to avoid computing xddt[a,b] several times.
            Integer temp = 0; 
            for (const BinWord &x : xddt_a_b ){
                if (cpp_scal_prod(u,x)^cpp_scal_prod(v,s[x])==0){
                    temp ++;
                }
                else {
                    temp --;
                }
            }
            row.push_back(temp);
        }
        result.push_back(row);
    }
    return result;
}