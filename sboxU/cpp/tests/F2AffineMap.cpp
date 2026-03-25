#include "../sboxU-cpp-module.hpp"

void print_cpp_F2AffineMap(const cpp_F2AffineMap &A){
    std::string res = "[";
    for (auto b : cpp_to_bin(A.cstte,A.get_output_length())){
        res = res + to_string(b) + " ";
    }
    res = res + "] \n + \n";
    for (auto v : A.get_image_vectors()){
        for (auto b : cpp_to_bin(v,A.get_input_length())){
        res = res + to_string(b) + " ";
    }
    res += '\n';
    }
    res+='\n';
    std::cout << res;
}

int main()
{   cpp_F2AffineMap A=cpp_F2AffineMap({1,2,4,8},4,4,0);
    print_cpp_F2AffineMap(A);
    cpp_F2AffineMap B=cpp_F2AffineMap({8,4,2,1},4,4,0);
    print_cpp_F2AffineMap(B);
    print_cpp_F2AffineMap(A+B);
    print_cpp_F2AffineMap(8+A+B);
    cpp_F2AffineMap C= 3 + A;
    for (auto y : C.get_cpp_S_box().get_lut()){
        std::cout << y << ", ";
    }
    std::cout << '\n' << C.get_cstte() <<'\n';
}