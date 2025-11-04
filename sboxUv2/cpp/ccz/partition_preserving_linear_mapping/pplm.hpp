//
// Created by Jules Baudrin on 19/06/2025.
// Modified by Jules Baudrin on 03/11/2025.
//

#ifndef PARTITION_PRESERVING_LINEAR_MAPPING_HPP
#define PARTITION_PRESERVING_LINEAR_MAPPING_HPP

#include <iostream>
#include <omp.h>

#include <bit>
#include <cmath>
#include <algorithm>
#include <vector>
#include <map>
#include <set>
#include <cassert>
#include <stdexcept>

#include "../../core/include.hpp"
#include "../../statistics/include.hpp"

#define BINWORD_SIZE 64
#define EMPTY UINT64_MAX
#define ONE ((BinWord) 1)

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
 *
 * Some methods are also described:
 * - in_which_bucket returns the key corresponding to the only bucket in which the parameter belongs to
 * - adapted_basis selects a basis of GF(2**n) suited to the search for partition-preserving linear mappings.
 */
template <class keyType>
class Partition {
private:
    uint64_t dimension; // Dimension of the space.
    map<keyType, vector<BinWord>> partition; // A partition of the space according to some key type.
    map<BinWord, keyType> reverse_dict; // The reverse dictionary.
public:

    vector<BinWord> basis; // A basis of the space.

    Partition(map<keyType, vector<BinWord>> _partition, Integer _dimension) : dimension(_dimension), partition(_partition)
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
                if(cpp_rank_of_vector_set(basis) != (Integer)basis.size())
                    basis.pop_back();
            }
        }
    }

    bool key_exists(const keyType &key) const{
        if(partition.count(key) == 0)
            return false;
        else
            return true;
    }

    bool is_compatible_with(const Partition<keyType> &other) const {
        for(const auto &[k, v]: partition) {
            if(other.key_exists(k)) {
                if(v.size() != other[k].size())
                    return false;
            }
            else
                return false;
        }
        return true;
    }

    vector<BinWord>& operator[] (keyType key) {return partition[key];}
    const vector<BinWord>& operator[](keyType key) const {return partition.at(key);}
};

bool cpp_verify_el_equivalence(cpp_S_box *sbox1,  cpp_S_box *sbox2, const cpp_BinLinearMap &solution);

class PPExitCondition
{
public:
    virtual bool is_valid(const cpp_BinLinearMap &solution) = 0;
};

class PPExitConditionTrivial : public virtual PPExitCondition
{
public:
  bool is_valid(const cpp_BinLinearMap &solution) {return true;}
};

class PPExitConditionEL : public virtual PPExitCondition
{
private:
    cpp_S_box *sbox1, *sbox2;

public:
    PPExitConditionEL(cpp_S_box *_sbox1, cpp_S_box *_sbox2) {
        sbox1 = _sbox1;
        sbox2 = _sbox2;
    }

    bool is_valid(const cpp_BinLinearMap &solution) {
        if(cpp_verify_el_equivalence(sbox1, sbox2, solution) == false)
            cout << "ICI PB " << endl;
        return cpp_verify_el_equivalence(sbox1, sbox2, solution);
    }
};


// ======================
class PPEarlyAbortCondition
{
public:
    virtual bool is_valid(const BinWord &x, const BinWord &y, const uint64_t &n) = 0;
};

class PPEarlyAbortConditionTrivial : public virtual PPEarlyAbortCondition
{
public:
    bool is_valid(const BinWord &x, const BinWord &y, const uint64_t &n) {return true;}
};

class PPEarlyAbortConditionEL : public virtual PPEarlyAbortCondition
{
public:
    bool is_valid(const BinWord &pre_image, const BinWord &image, const uint64_t &n) {
        BinWord mask = 0;
        for (int i = 0; i < n; ++i)
            mask |= ((BinWord)1) << i;

        BinWord pre_image_y = (pre_image >> n) & mask;
        BinWord image_y = (image >> n) & mask;

        if(pre_image_y)
            return true;
        else{
            if(!image_y)
                return true;
            else {
                return false;
            }
        }
    }
};


/*
 * ===========================================
 * The PartitionPreservingLinearMappings class
 * ===========================================
 * It describes a search strategy for linear bijections which send an input partition to an output one. It is essentially made of :
 * - two pointers to the input/output partitions.
 * - the Boolean variable single_non_trivial_answer indicates whether a full search is perfomed or not.
 * - the integer number_of_threads.
 * - a vector of solutions found_mappings.
 *
 * Some methods are also described:
 * - recursive_search handles the search. It follows a DFS strategy.
 * - search launches the recursive search and handles the parallelization.
 */
template <class keyType>
class PartitionPreservingLinearMappings {
private:
    bool single_non_trivial_answer;
    unsigned int number_of_threads;
    uint64_t dimension;
    Partition<keyType> *partition_in;
    Partition<keyType> *partition_out;
    PPExitCondition *exit_condition;
    PPEarlyAbortCondition *early_abort;

public:
    vector<cpp_BinLinearMap> found_mappings;

    PartitionPreservingLinearMappings(Partition<keyType> *_partition_in, Partition<keyType> *_partition_out, PPExitCondition *_exit_condition, PPEarlyAbortCondition *_early_abort, bool _single_non_trivial_answer, unsigned int _number_of_threads) :
        single_non_trivial_answer(_single_non_trivial_answer),
        number_of_threads(_number_of_threads),
        partition_in(_partition_in),
        partition_out(_partition_out),
        exit_condition(_exit_condition),
        early_abort(_early_abort) {
        assert(_partition_in->dim() == _partition_out->dim());
        dimension = _partition_in->dim();
    }

    vector<cpp_BinLinearMap> recursive_search(vector<BinWord> basis_out) {
        unsigned int current_index = basis_out.size();

        // A leaf is reached
        if(current_index == dimension) {
            // Computes the corresponding BinLinearMap
            vector<BinWord> image_canonical_basis(dimension);
            BinWord x, y;
            for(BinWord i = 0; i < (ONE << dimension); i++) {
                x = cpp_linear_combination(partition_in->basis, i);
                if(popcount(x) == 1) { // if x is a vector of the canonical basis
                    y = cpp_linear_combination(basis_out, i);
                    image_canonical_basis[countr_zero(x)] = y; // countr_zero provides the index of the first set bit.
                }
            }
            cpp_BinLinearMap solution(image_canonical_basis); 
            solution = solution.transpose().inverse();// CAREFUL HERE


            // If the exit condition is not valid or if a single solution is wanted but this one is the identity, returns nothing
            if(!(*exit_condition).is_valid(solution) || (single_non_trivial_answer && basis_out == partition_in->basis)) 
                return {};
            else // otherwise, return the mapping
                return {solution};
        }

        vector<cpp_BinLinearMap> solutions; // The list of the found solutions
        BinWord x, y, v, w, image, pre_image, i;
        bool partition_preserving = true;

        x = partition_in->basis[current_index]; // The preimage of the currently guessed image.
        const keyType bucket_key = partition_in->in_which_bucket(x); // The bucket to inspect
        const vector<BinWord> bucket = (*partition_out)[bucket_key];

        for(unsigned int k = 0; k < bucket.size() && !(single_non_trivial_answer && solutions.size()); k++) { // For each possible image
            y = bucket[k];

            basis_out.push_back(y);

            // Early abort if the dimension does not increase (bijective mappings)
            if(cpp_rank_of_vector_set(basis_out) == (Integer) (current_index + 1)) {

                // The partition-preserving condition : each image must belong to the out_bucket associated to in_bucket to which the preimage belongs.
                // The condition is only tested on (preimage, image) pairs that have not been tested before
                for(i = 0; i < (1 << (current_index)) && partition_preserving; i++) {
                    v = cpp_linear_combination(partition_in->basis, i); // The i-th vector of Span(basis_in)
                    w = cpp_linear_combination(basis_out, i); // The i-th vector of Span(basis_out)
                    pre_image = x ^ v;
                    image = y ^ w;

                    if(partition_out->in_which_bucket(image) != partition_in->in_which_bucket(pre_image) || !(*early_abort).is_valid(pre_image, image, dimension >> 1))
                        partition_preserving = false;
                }

                // Recursive call if no contradiction for now
                if(partition_preserving) {
                    const vector<cpp_BinLinearMap> recursive_solutions = recursive_search(basis_out);
                    solutions.insert(solutions.end(), recursive_solutions.begin(), recursive_solutions.end());  // Adds the found mappings to the list
                }
            }
            basis_out.pop_back();
            partition_preserving = true;
        }
        return solutions;
    }

    void search() {
        // Early abort if the partitions are not compatible
        if(!partition_in->is_compatible_with(*partition_out))
            return;

        if(number_of_threads == 1) {
            found_mappings = recursive_search({});
        }
        else {
            keyType bucket_key;
            bucket_key = partition_in->in_which_bucket(partition_in->basis[0]);
            const vector<BinWord> bucket = (*partition_out)[bucket_key];

            bool found = false;
            vector<vector<cpp_BinLinearMap>> thread_solutions(number_of_threads);

            omp_set_num_threads(number_of_threads);
#pragma omp parallel for default(none) shared(bucket, found, thread_solutions, cout)
            for(unsigned int i = 0; i < bucket.size(); i++) {
                if(!(single_non_trivial_answer && found)) {
                    vector<cpp_BinLinearMap> cur_solutions;
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

vector<cpp_BinLinearMap> cpp_equivalences_from_lat(
    cpp_S_box sbox1,
    cpp_S_box sbox2,
    const bool &single_non_trivial_answer,
    const unsigned int &number_of_threads,
    const string equivalence_type);

#endif // PARTITION_PRESERVING_LINEAR_MAPPING_HPP
