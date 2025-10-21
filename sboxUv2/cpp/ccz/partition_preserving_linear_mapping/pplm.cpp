//
// Created by Jules Baudrin on 19/06/2025.
// Modified by Jules Baudrin on 21/10/2025.
//
#include "pplm.hpp"

vector<vector<cpp_BinLinearMap>> cpp_equivalences_from_lat(
    const vector<vector<Integer>> &lat1,
    const vector<vector<Integer>> &lat2,
    const bool &single_non_trivial_answer,
    const unsigned int &number_of_threads,
    const string equivalence_type)
{
	Integer n = log2(lat1.size());
	BinWord size = (ONE << n);

	map <pair<Integer, string>, vector<BinWord>> dict1;
	map <pair<Integer, string>, vector<BinWord>> dict2;
	for(BinWord x = 0; x < size; x++) {
		for(BinWord y = 0; y < size; y++) {
			if(equivalence_type == "ccz-linear" || (x != 0 && y != 0)) { // Either we search for ccz equivalences and all (x, y) are mapped in "mixed" buckets
				dict1[{lat1[x][y], "mixed"}].push_back((y << n) ^ x);
				dict2[{lat2[x][y], "mixed"}].push_back((y << n) ^ x);
			}
			else { // Or we search for EL or L equivalences
				if(x == 0) {// Otherwise, split the buckets so that (x, 0) is necessarily mapped onto (x', 0) (and the same for (0, y) -> (0, y') in the case of linear equivalence).
					dict1[{lat1[x][y], "x=0"}].push_back((y << n) ^ x);
					dict2[{lat2[x][y], "x=0"}].push_back((y << n) ^ x);
				}
				else { // y == 0
					if(equivalence_type == "linear") {
						dict1[{lat1[x][y], "y=0"}].push_back((y << n) ^ x);
						dict2[{lat2[x][y], "y=0"}].push_back((y << n) ^ x);
					}
					else {
						dict1[{lat1[x][y], "mixed"}].push_back((y << n) ^ x);
						dict2[{lat2[x][y], "mixed"}].push_back((y << n) ^ x);
					}
				}
			}
		}
	}

	Partition<pair<Integer, string>> partition1(dict1, 2 * n);
	Partition<pair<Integer, string>> partition2(dict2, 2 * n);

	PartitionPreservingLinearMappings<pair<Integer, string>> pplm(&partition1, &partition2, single_non_trivial_answer, number_of_threads);
	pplm.search();

	vector<vector<cpp_BinLinearMap>> block_decompositions;
	for(const auto &lut: pplm.found_mappings) {
		// A D
		// C B

		vector<vector<BinWord>> image_vectors(4); // A B C D
		BinWord v, w;
		for (int i = 0; i < n; i++) {
			v = lut[1 << i];
			w = lut[1 << (n + i)];
			image_vectors[0].push_back((v << n) >> n); //A
			image_vectors[2].push_back(v >> n); // C

			image_vectors[3].push_back((w << n) >> n); // D
			image_vectors[1].push_back(w >> n); // B
		}
		vector<cpp_BinLinearMap> block_dec;
		for (int i = 0; i < image_vectors.size(); i++) {
			block_dec.push_back(cpp_BinLinearMap(image_vectors[i], n, n));
		}
		block_decompositions.push_back(block_dec);
	}
	return block_decompositions;
}
