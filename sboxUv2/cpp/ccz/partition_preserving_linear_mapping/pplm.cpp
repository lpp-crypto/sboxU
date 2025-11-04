//
// Created by Jules Baudrin on 19/06/2025.
// Modified by Jules Baudrin on 03/11/2025.
//
#include "pplm.hpp"


bool cpp_verify_el_equivalence(cpp_S_box *sbox1,  cpp_S_box *sbox2, const cpp_BinLinearMap &solution) {
    vector<cpp_BinLinearMap> ABCD = cpp_ccz_block_decomposition(solution);
 	cpp_S_box A = ABCD[0].get_cpp_S_box();
    cpp_S_box B = ABCD[1].get_cpp_S_box();
    cpp_S_box C = ABCD[2].get_cpp_S_box();
    cpp_S_box D = ABCD[3].get_cpp_S_box();
    
    for (int i = 0; i < sbox1->size(); ++i) {
    	if(D[i])
    		return false;
    }
    auto sbox1_prime = (B * (*sbox1)) + C;
    auto sbox2_prime = (*sbox2) * A;
    if (sbox1_prime != sbox2_prime)
    	cout << "EL equivalence test failed" << endl;
    return sbox1_prime == sbox2_prime;
}


vector<cpp_BinLinearMap> cpp_equivalences_from_lat(
    cpp_S_box sbox1,
    cpp_S_box sbox2,
    const bool &single_non_trivial_answer,
    const unsigned int &number_of_threads,
    const string equivalence_type)
{
	Integer n = log2(sbox1.size());
	BinWord size = (ONE << n);
	vector<vector<Integer>> lat1 = cpp_lat(sbox1);
	vector<vector<Integer>> lat2 = cpp_lat(sbox2);

	map<pair<Integer, string>, vector<BinWord>> dict1;
	map<pair<Integer, string>, vector<BinWord>> dict2;
	for(BinWord x = 0; x < size; x++) {
		for(BinWord y = 0; y < size; y++) {
			if(equivalence_type == "ccz-linear" || equivalence_type == "extended-linear") { // Either we search for ccz equivalences and all (x, y) are mapped in "mixed" buckets
				dict1[{lat1[x][y], "mixed"}].push_back((y << n) ^ x);
				dict2[{lat2[x][y], "mixed"}].push_back((y << n) ^ x);
			}
			if(equivalence_type == "linear") {
				if(x == 0) { 
					dict1[{lat1[x][y], "x=0"}].push_back((y << n) ^ x);
					dict2[{lat2[x][y], "x=0"}].push_back((y << n) ^ x);
				}
				else if(y == 0) {
					dict1[{lat1[x][y], "y=0"}].push_back((y << n) ^ x);
					dict2[{lat2[x][y], "y=0"}].push_back((y << n) ^ x);
				}
				else {
					dict1[{lat1[x][y], "mixed"}].push_back((y << n) ^ x);
					dict2[{lat2[x][y], "mixed"}].push_back((y << n) ^ x);
				}
			}
			// if(equivalence_type == "extended-linear") {
			// 	if(y == 0) { 
			// 		dict1[{lat1[x][y], "y=0"}].push_back((y << n) ^ x);
			// 		dict2[{lat2[x][y], "y=0"}].push_back((y << n) ^ x);
			// 	}
			// 	else {
			// 		dict1[{lat1[x][y], "mixed"}].push_back((y << n) ^ x);
			// 		dict2[{lat2[x][y], "mixed"}].push_back((y << n) ^ x);
			// 	}
			// }
		}
	}    

	Partition<pair<Integer, string>> partition1(dict1, 2 * n);
	Partition<pair<Integer, string>> partition2(dict2, 2 * n);

	PPExitCondition *exit_condition; 
	PPEarlyAbortCondition *early_abort_condition;
	if(equivalence_type == "extended-linear") {
		exit_condition = new PPExitConditionEL(&sbox1, &sbox2);
		early_abort_condition = new PPEarlyAbortConditionEL();
	}
	else {
		exit_condition = new PPExitConditionTrivial();
		early_abort_condition = new PPEarlyAbortConditionTrivial();
	}
	
		
	PartitionPreservingLinearMappings<pair<Integer, string>> pplm = PartitionPreservingLinearMappings<pair<Integer, string>>(&partition1, &partition2, exit_condition, early_abort_condition, single_non_trivial_answer, number_of_threads);
	pplm.search();
	return pplm.found_mappings;
}
