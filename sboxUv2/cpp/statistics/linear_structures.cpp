#include "linear_structures.hpp"


std::vector<std::vector<Integer> > cpp_linear_structures(const cpp_S_box & f){
    std::vector<std::vector<Integer>> result;
    std::vector<Integer> l0;
    result.push_back(l0);
    std::vector<Integer> l1;
    result.push_back(l1);
    for (Integer a = 1; a < f.input_space_size(); a++) {
        Integer offset=f[a]^f[0];
        bool valid=true;
        for (Integer x = 1; x <f.input_space_size();x++){
            if (f[x^a]^f[x] != offset) {
                valid=false;
                break;
            }

        }
        if (valid) {
            result[offset].push_back(a);
        }
    }
    return result;
}

std::map<Integer,std::vector<std::vector<Integer> >> cpp_linear_structures_vectorial(const cpp_S_box &s){
    std::map<Integer,std::vector<std::vector<Integer> >> result;
    for (Integer c = 1; c < s.input_space_size(); c++){
        cpp_S_box f=s.component(c);
        std::vector<std::vector<Integer> > l = cpp_linear_structures(f);
        if (l[0].size()>0 or l[1].size()>0){
            result.insert({c,l});
        }
    
    }
    return result;
}

cpp_Spectrum cpp_linear_structures_vectorial_spectrum(const cpp_S_box &s){
    cpp_Spectrum result;
    std::map<Integer,std::vector<std::vector<Integer> >> l = cpp_linear_structures_vectorial(s);
    for (auto c : l){
    result.incr_by_amount(c.first,c.second[0].size()+c.second[1].size());
    }
    return result;
}

