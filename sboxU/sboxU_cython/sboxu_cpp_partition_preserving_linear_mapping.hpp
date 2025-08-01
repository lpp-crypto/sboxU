// Created on 2025-06-19 by baudrin-j.
// Modified on 2025-07-27 by baudrin-j.

#ifndef PARTITION_PRESERVING_LINEAR_MAPPING_HPP
#define PARTITION_PRESERVING_LINEAR_MAPPING_HPP

#include <iostream>
#include <omp.h>

#include <cmath>
#include <algorithm>
#include <vector>
#include <map>
#include <set>
#include <cassert>
#include <stdexcept>

#include "sboxu_cpp_utils.hpp"

#define BINWORD_SIZE 64

using namespace std;

/*
 * =======================
 * The Partition class
 * =======================
 * It describes a partition of GF(2)^n. It is made of :
 * - a dictionary whose key can be anything and whose value is a subset of GF(2)^n (represented as a vector of BinWord).
 *      The subsets form a partition of GF(2)^n.
 * - a "reverse dictionary" which maps any value of GF(2)^n to the corresponding key.
 * - an int corresponding to the dimension of GF(2)^n, i.e equal to n.
 * - a basis of GF(2)^n. (or of GF(2)^(n/2) in some cases).
 * - (in some cases another basis of GF(2)^(n/2))
 *
 * Some methods are also described:
 * - in_which_bucket returns the key corresponding to the only bucket in which the parameter belongs to
 * - adapted_basis (and diagonal_adapted_basis) selects a basis of GF(2**n) (two bases of GF(2**(n/2))) suited to the search for (diagonal) partition-preserving linear mappings.
 */
template <class keyType>
class Partition {
private:
	uint64_t dimension; // Dimension of the space.
	map<keyType, vector<BinWord>> partition; // A partition of the space according to some key type.
	map<BinWord, keyType> reverse_dict; // The reverse dictionary.
public:

	vector<BinWord> basis; // A basis of the space.
	vector<BinWord> basis_B; // only used in case of diagonal mappings

	Partition(map<keyType, vector<BinWord>> _partition, Integer _dimension, bool diagonal_mappings) : dimension(_dimension), partition(_partition)
	{
		set<BinWord> all_values;
		for(const auto &[k, v]: partition) {
			Integer previous_size = all_values.size();
			all_values.insert(v.begin(), v.end());
			assert(previous_size + v.size() == all_values.size()); // All sets must be disjoint.
			for(const auto &x: v) {
				assert(x < (ONE << dimension)); // Each element belongs to [0, 2^n-1]
				reverse_dict[x] = k;
			}
		}
		assert(all_values.size() == ONE << dimension); // all_values = [0, 2^n-1]
		if(diagonal_mappings)
			diagonal_adapted_basis();
		else
			adapted_basis();
	};

	uint64_t dim() {return dimension;}

	// Return the key corresponding to the only bucket to which elt belongs
	keyType in_which_bucket(const BinWord &elt) {
		return reverse_dict[elt];
	}

	// Return a sorted list of pair (length, key) which corresponds to the ordering of the partition {k:v} by increasing length of its values (v).
	vector<pair<int, keyType>> _increasing_order() {
		vector<pair<int, keyType>> size_key_table;

		for (const auto &v : partition)
			size_key_table.push_back({v.second.size(), v.first}); // (size, key) pair
		sort(size_key_table.begin(), size_key_table.end());

		return size_key_table;
	}

	// Select a better basis which minimizes the number of choices when looking for partition-preserving linear mappings.
	void adapted_basis() {
		basis.clear();
		// Enumerate the partition {k : v} by increasing length of v
		for(const auto &[len, key]: _increasing_order()) {
			if(basis.size() == dimension)
				break;
			// Enumerate all elements of v
			for(const auto &x: partition[key]) {
				// Add x to the basis only if the rank increases
				basis.push_back(x);
				if(rank_of_vector_set_cpp(basis) != (Integer)basis.size())
					basis.pop_back();
			}
		}
	}

	// Select an *independent family of size n/2* which minimizes the number of choices when looking for partition-preserving linear mappings that are DIAGONAL.
	// The family {f_i} has the property that Span(f_i >> n/2, 0 <= i <= dim/2) = GF(2**n/2)  and Span((f_i << n/2) >> n/2, 0 <= i <= dim/2)= GF(2**n/2).
	void diagonal_adapted_basis() {
		assert(dimension % 2 == 0);
		uint64_t half_dim = dimension / 2;

		basis.clear(); // A basis of GF(2**(n/2))
		basis_B.clear(); // Another basis of GF(2**(n/2))
		// Enumerate the partition {k : v} by increasing length of v
		for(const auto &[len, key]: _increasing_order()) {
			if(basis.size() == half_dim)
				break;
			// Enumerate all elements of v
			for(const auto &x: partition[key]) {
				// Add x to the basis only if the rank increases
				basis.push_back((x << (BINWORD_SIZE - half_dim)) >> (BINWORD_SIZE - half_dim));
				basis_B.push_back((x << (BINWORD_SIZE - dimension)) >> (BINWORD_SIZE - half_dim));
				if(rank_of_vector_set_cpp(basis) != (Integer)basis.size() || rank_of_vector_set_cpp(basis_B) != (Integer)basis_B.size()) {
					basis.pop_back();
					basis_B.pop_back();
				}
			}
		}
	}

	vector<BinWord>& operator[] (keyType key) {return partition[key];}
	const vector<BinWord>& operator[](keyType key) const {return partition.at(key);}
};


/*
 * ===========================================
 * The PartitionPreservingLinearMappings class
 * ===========================================
 * It describes a search strategy for linear bijections which send an input partition to an output one. It is essentially made of :
 * - two pointers to the input/output partitions.
 * - the Boolean variable single_non_trivial_answer indicates whether a full search is perfomed or not.
 * - the Boolean variable diagonal_mappings indicates whether the linear bijections must have the form diag(A, B) or not.
 * - the integer number_of_threads.
 * - a vector of solutions found_mappings.
 *
 * Some methods are also described:
 * - recursive_search and diagonal_recursive_search which handles the search. They follow the same DFS strategy.
 * diagonal_recursive_search only considers mappings of the form diag(A, B).
 * - recursive_partition_preserving_linear_mappings is the recursive search for partition-preserving linear mappings.
 */
template <class keyType>
class PartitionPreservingLinearMappings {
private:
	bool single_non_trivial_answer;
	bool diagonal_mappings;
	unsigned int number_of_threads;
	uint64_t dimension;
	Partition<keyType> *partition_in;
	Partition<keyType> *partition_out;

public:
	vector<vector<BinWord>> found_mappings;

	PartitionPreservingLinearMappings(Partition<keyType> *_partition_in, Partition<keyType> *_partition_out, bool _single_non_trivial_answer, bool _diagonal_mappings, unsigned int _number_of_threads) :
			single_non_trivial_answer(_single_non_trivial_answer),
			diagonal_mappings(_diagonal_mappings),
			number_of_threads(_number_of_threads),
			partition_in(_partition_in),
			partition_out(_partition_out) {
			assert(_partition_in->dim() == _partition_out->dim());

			if(!diagonal_mappings)
				dimension = _partition_in->dim();
			else {
				assert(_partition_in->dim() % 2 == 0);
				dimension = _partition_in->dim() / 2;
			}
		}

	vector<vector<BinWord>> recursive_search(vector<BinWord> basis_out) {
		unsigned int current_index = basis_out.size();
		// End condition, return the found solution
		if(current_index == dimension) {
			if(single_non_trivial_answer && basis_out == partition_in->basis) // If a single solution is wanted but this one is the identity, returns nothing
				return {};
			else { // In any other case, return the LUT
				vector<BinWord> lut(ONE << dimension);
				BinWord x, y;
				for(BinWord i = 0; i < (ONE << dimension); i++) {
					x = matrix_vector_multiplication(i, partition_in->basis);
					y = matrix_vector_multiplication(i, basis_out);
					lut[x] = y;
				}
				return {lut};
			}
		}

		vector<vector<BinWord>> solutions; // The list of the found solutions
		BinWord x, y, v, w, image, pre_image, i;
		bool partition_preserving = true;

		x = partition_in->basis[current_index]; // The preimage of the currently guessed image.
		const keyType bucket_key = partition_in->in_which_bucket(x); // The bucket to inspect
		const vector<BinWord> bucket = (*partition_out)[bucket_key];

		for(unsigned int k = 0; k < bucket.size() && !(single_non_trivial_answer && solutions.size()); k++) { // For each possible image
			y = bucket[k];

			basis_out.push_back(y);

			// Early abort if the dimension does not increase (bijective mappings)
			if(rank_of_vector_set_cpp(basis_out) == (Integer) (current_index + 1)) {

				// The partition-preserving condition : each image must belong to the out_bucket associated to in_bucket to which the preimage belongs.
				// The condition is only tested on (preimage, image) pairs that have not been tested before
				for(i = 0; i < (1 << (current_index)) && partition_preserving; i++) {
					v = matrix_vector_multiplication(i, partition_in->basis); // The i-th vector of Span(basis_in)
					w = matrix_vector_multiplication(i, basis_out); // The i-th vector of Span(basis_out)
					pre_image = x ^ v;
					image = y ^ w;

					if(partition_out->in_which_bucket(image) != partition_in->in_which_bucket(pre_image))
						partition_preserving = false;
				}

				// Recursive call if no contradiction for now
				if(partition_preserving) {
					const vector<vector<BinWord>> recursive_solutions = recursive_search(basis_out);
					solutions.insert(solutions.end(), recursive_solutions.begin(), recursive_solutions.end());  // Adds the found mappings to the list
				}
			}
			basis_out.pop_back();
			partition_preserving = true;
		}
		return solutions;
	}

	vector<vector<BinWord>> diagonal_recursive_search(vector<BinWord> basis_out_A, vector<BinWord> basis_out_B) {
		unsigned int current_index = basis_out_A.size();

		// End condition, return the found solution
		if(current_index == dimension) {
			if(single_non_trivial_answer && basis_out_A == partition_in->basis && basis_out_B == partition_in->basis_B) // If a single solution is wanted but this one is the identity, returns nothing
				return {};
			else { // In any other case, return the LUTs of A and B (the one of B is appended to the one of A).
				vector<BinWord> luts_AB(ONE << (dimension + 1));
				BinWord x, y;
				for(BinWord i = 0; i < (ONE << dimension); i++) {
					x = matrix_vector_multiplication(i, partition_in->basis);
					y = matrix_vector_multiplication(i, basis_out_A);
					luts_AB[x] = y;
					x = matrix_vector_multiplication(i, partition_in->basis_B);
					y = matrix_vector_multiplication(i, basis_out_B);
					luts_AB[(ONE << dimension) + x] = y;
				}
				return {luts_AB};
			}
		}

		vector<vector<BinWord>> solutions; // The list of the found solutions
		BinWord x_A, x_B, y_A, y_B, image, pre_image, i, j;
		bool partition_preserving = true;

		x_A = partition_in->basis[current_index]; // The preimage of the currently guessed image.
		x_B = partition_in->basis_B[current_index];
		const keyType bucket_key = partition_in->in_which_bucket((x_B << dimension) | x_A); // The bucket to inspect
		const vector<BinWord> bucket = (*partition_out)[bucket_key];

		for(unsigned int k = 0; k < bucket.size() && !(single_non_trivial_answer && found_mappings.size()); k++) { // For each possible image
			y_A = (bucket[k] << (BINWORD_SIZE - dimension)) >> (BINWORD_SIZE - dimension);
			y_B = (bucket[k] << (BINWORD_SIZE - dimension * 2)) >> (BINWORD_SIZE - dimension);

			basis_out_A.push_back(y_A);
			basis_out_B.push_back(y_B);

			// Early abort if one of the dimension does not increase
			if(rank_of_vector_set_cpp(basis_out_A) == (Integer) (current_index + 1) &&
			   rank_of_vector_set_cpp(basis_out_B) == (Integer) (current_index + 1)) {

				// The partition-preserving condition : each image must belong to the out_bucket associated to in_bucket to which the preimage belongs.
				// The condition is only tested on (preimage, image) pairs that have not been tested before
				for(i = 0; i < (1 << (current_index + 1)) && partition_preserving; i++) {
					x_A = matrix_vector_multiplication(i, partition_in->basis); // The i-th vector of <b for b in basis_in>
					y_A = matrix_vector_multiplication(i, basis_out_A); // The i-th vector of <b for b in basis_ou_At>
					for(j = 0; j < (1 << (current_index + 1)) && partition_preserving; j++) {
						x_B = matrix_vector_multiplication(j, partition_in->basis_B); // The j-th vector of <b for b in basis_B>
						y_B = matrix_vector_multiplication(j, basis_out_B); // The j-th vector of <b for b in basis_out_B>

						if(i > (1 << current_index) || j > (1 << current_index)) {
							pre_image = (x_B << dimension) | x_A;
							image = (y_B << dimension) | y_A;
							if(partition_out->in_which_bucket(image) != partition_in->in_which_bucket(pre_image))
								partition_preserving = false;
						}
					}
				}

				// Recursive call if no contradiction for now
				if(partition_preserving) {
					const vector<vector<BinWord>> recursive_solutions = diagonal_recursive_search(basis_out_A, basis_out_B);
					solutions.insert(solutions.end(), recursive_solutions.begin(), recursive_solutions.end());  // Adds the found mappings to the list
				}
			}
			basis_out_A.pop_back();
			basis_out_B.pop_back();
			partition_preserving = true;
		}
		return solutions;
	}

	void search() {
		if(number_of_threads == 1) {
			if(diagonal_mappings)
				found_mappings = diagonal_recursive_search({}, {});
			else
				found_mappings = recursive_search({});
		}
		else {
			keyType bucket_key;
			if(diagonal_mappings)
				bucket_key = partition_in->in_which_bucket((partition_in->basis_B[0] << dimension) | partition_in->basis[0]);
			else
				bucket_key = partition_in->in_which_bucket(partition_in->basis[0]);
			const vector<BinWord> bucket = (*partition_out)[bucket_key];

			bool found = false;
			vector<vector<vector<BinWord>>> thread_solutions(number_of_threads);

			omp_set_num_threads(number_of_threads);
		#pragma omp parallel for default(none) shared(bucket, found, thread_solutions, cout)
			for(unsigned int i = 0; i < bucket.size(); i++) {
				if(!(single_non_trivial_answer && found)) {
					vector <vector<BinWord>> cur_solutions;
					if(diagonal_mappings) {
						const BinWord y_A = (bucket[i] << (BINWORD_SIZE - dimension)) >> (BINWORD_SIZE - dimension);
						const BinWord y_B = (bucket[i] << (BINWORD_SIZE - dimension * 2)) >> (BINWORD_SIZE - dimension);
						cur_solutions = diagonal_recursive_search({y_A}, {y_B});
					} else
						cur_solutions = recursive_search({bucket[i]});
					if(single_non_trivial_answer && cur_solutions.size()) {
						#pragma omp critical
						found = true;
					}
					thread_solutions[omp_get_thread_num()].insert(thread_solutions[omp_get_thread_num()].end(), cur_solutions.begin(), cur_solutions.end());
				}
			}

			for(auto & v : thread_solutions){
				found_mappings.insert(found_mappings.end(), v.begin(), v.end());
			}
		}
	}
};

vector<pair<vector<BinWord>, vector<BinWord>>> linear_automorphism_search(const vector<vector<Integer>> &lat,  const bool &single_non_trivial_answer, const bool &alternative_partition, const unsigned int &number_of_threads);

pair<vector<BinWord>, vector<BinWord>> cpp_is_linearly_self_equivalent_from_lat(const vector<vector<Integer>> &lat, const string algo, const unsigned int &number_of_threads);
vector<pair<vector<BinWord>, vector<BinWord>>> cpp_linear_automorphisms_from_lat(const vector<vector<Integer>> &lat, const string algo, const unsigned int &number_of_threads);

#endif // PARTITION_PRESERVING_LINEAR_MAPPING_HPP