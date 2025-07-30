//
// Created on 2025-06-19 by baudrin-j.
// Modified on 2025-07-27 by baudrin-j.
//
#include "sboxu_cpp_partition_preserving_linear_mapping.hpp"

// The core search function
//  1) single_non_trivial_answer
//      - if true, the search is stopped whenever non-trivial automorphism is found (if found).
//		- otherwise, a full search is performed.
//  2) diagonal_mappings
//      - if true, the automorphisms are built under the form diag(A, B). Can be used both alternative_partition=true or false.
//      - otherwise, generic 2n-bit-mappings are searched. The combination (diagonal_mappings = false, alternative_partition=true) still
//        returns automorphisms of the form diag(A, B). In that case, the constraint comes from the partition.
//	3) alternative_partition
//      - if true, the standard partition is used : {lat[x][y], x, y \in GF(2**n)} is split by values.
//      - otherwise, {lat[x][y], x, y \in GF(2**n)} is split by values. Furthermore, the buckets are split
//        so that (x, 0) is necessarily mapped onto (x', 0) (and the same for (0, y) -> (0, y') ).
vector<pair<vector<BinWord>, vector<BinWord>>> linear_automorphism_search(const vector<vector<Integer>> &lat, const bool &single_non_trivial_answer, const bool &diagonal_mappings, const bool &alternative_partition, const unsigned int &number_of_threads) {
	Integer n = log2(lat.size());
	BinWord size = (ONE << n);

	map <pair<Integer, string>, vector<BinWord>> dict;
	for(BinWord x = 0; x < size; x++) {
		for(BinWord y = 0; y < size; y++) {
			if(!alternative_partition || (x != 0 && y != 0)) // if !alternative_partition, only take into account the value lat[x][y] ('mixed' is appended for all x, y)
				dict[{lat[x][y], "mixed"}].push_back((y << n) ^ x);
			else if(x == 0) // Otherwise, split the buckets so that (x, 0) is necessarily mapped onto (x', 0) (and the same for (0, y) -> (0, y') ).
				dict[{lat[x][y], "y only"}].push_back((y << n) ^ x);
			else // y == 0
				dict[{lat[x][y], "x only"}].push_back((y << n) ^ x);
		}
	}
	Partition<pair<Integer, string>> partition(dict, 2 * n, diagonal_mappings);
	PartitionPreservingLinearMappings<pair<Integer, string>> pplm(&partition, &partition, single_non_trivial_answer, diagonal_mappings, number_of_threads);
	pplm.search();

	vector<pair<vector<BinWord>, vector<BinWord>>> splitted_luts;
	for(const auto &lut: pplm.found_mappings) {
		vector<BinWord> lutA, lutB;
		for(uint i = 0; i < (ONE << n); i++) {
			lutA.push_back(lut[i]);
			if(diagonal_mappings) // If diagonal_mappings, lut contains the concatenation of lutA and lutB
				lutB.push_back(lut[(ONE << n) + i]);
			else // Otherwise, lutA and lutB must be extracted from the lut of diag(A, B)
				lutB.push_back(lut[i << n] >> n);
		}
		splitted_luts.push_back({lutA, lutB});
	}
	return splitted_luts;
}

// The three different search strategies
vector<pair<vector<BinWord>,vector<BinWord>>> linear_automorphisms(const vector<vector<Integer>> &lat, const bool &single_non_trivial_answer, const string &algo, const unsigned int &number_of_threads) {
	vector<pair<vector<BinWord>, vector<BinWord>>> luts;
	if(algo == "alt_partition_diag_mappings")
		luts = linear_automorphism_search(lat, single_non_trivial_answer, true, true, number_of_threads);
	else if(algo == "alt_partition")
		luts = linear_automorphism_search(lat, single_non_trivial_answer, false,true, number_of_threads);
	else if(algo == "std_partition_diag_mappings")
		luts = linear_automorphism_search(lat, single_non_trivial_answer, true, false, number_of_threads);

	return luts;
}

// The full search
vector<pair<vector<BinWord>, vector<BinWord>>> cpp_linear_automorphisms_from_lat(const vector<vector<Integer>> &lat, const string algo, const unsigned int &number_of_threads) {
	return linear_automorphisms(lat, false, algo, number_of_threads);
}

// The early-abort search
pair<vector<BinWord>, vector<BinWord>> cpp_is_linearly_self_equivalent_from_lat(const vector<vector<Integer>> &lat, const string algo, const unsigned int &number_of_threads) {
	auto luts = linear_automorphisms(lat, true, algo, number_of_threads);
	if(luts.size())
		return luts[0];
	else
		return {{}, {}};
}