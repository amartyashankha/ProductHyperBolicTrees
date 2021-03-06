#ifndef PRODUCT_TREES_HELPER_H
#define PRODUCT_TREES_HELPER_H

#include<iostream>
#include<vector>
#include<numeric>
#include<unordered_map>
#include<unordered_set>

#include<cassert>
#include<cmath>

#include <sparsehash/sparse_hash_map>
#include <sparsehash/dense_hash_map>
#include <sparsehash/sparse_hash_set>
#include <sparsehash/dense_hash_set>

using namespace std;

typedef unsigned int LONG_UINT;
//typedef google::dense_hash_map<LONG_UINT, unsigned int> UNORDERED_MAP;
typedef unordered_map<LONG_UINT, unsigned int> UNORDERED_MAP;

tuple<vector<LONG_UINT>, UNORDERED_MAP> get_subtree_labels(
        LONG_UINT *labels, LONG_UINT num_labels, unsigned int tree_depth, unsigned int subtree_depth);

tuple< LONG_UINT *, unsigned int * > create_subtree_partitions(
        LONG_UINT *labels, LONG_UINT num_labels,
        unsigned int tree_depth, unsigned int subtree_depth,
        vector<LONG_UINT> subtree_labels, UNORDERED_MAP map_subtree_label_to_index);

vector< pair<LONG_UINT, LONG_UINT> > get_all_edges(
        LONG_UINT *labels, LONG_UINT num_labels, unsigned int threshold, unsigned int tree_depth, unsigned int subtree_depth,
        vector<LONG_UINT> subtree_labels, UNORDERED_MAP map_subtree_label_to_index,
        LONG_UINT *subtree_label_to_label_index_list, unsigned int *subtree_member_list_indices);

#endif
