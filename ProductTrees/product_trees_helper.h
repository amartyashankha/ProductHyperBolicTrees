#ifndef PRODUCT_TREES_HELPER_H
#define PRODUCT_TREES_HELPER_H

#include<iostream>
#include<vector>
#include<numeric>
#include<unordered_map>

#include<cassert>
#include<cmath>

using namespace std;

typedef unsigned int LONG_UINT;
typedef unordered_map<LONG_UINT, unsigned int> UNORDERED_MAP;

tuple< vector<LONG_UINT>, unordered_map<LONG_UINT, unsigned int> > get_subtree_labels(
        LONG_UINT *labels, LONG_UINT num_labels, unsigned int tree_depth, unsigned int subtree_depth);

tuple< LONG_UINT *, unsigned int * > create_subtree_partitions(
        LONG_UINT *labels, LONG_UINT num_labels,
        unsigned int tree_depth, unsigned int subtree_depth,
        vector<LONG_UINT> subtree_labels, unordered_map<LONG_UINT, unsigned int> map_subtree_label_to_index);

vector< pair<LONG_UINT, LONG_UINT> > get_all_edges(
        LONG_UINT *labels, unsigned int threshold, unsigned int tree_depth, unsigned int subtree_depth,
        vector<LONG_UINT> subtree_labels, unordered_map<LONG_UINT, unsigned int> map_subtree_label_to_index,
        LONG_UINT *subtree_label_to_label_index_list, unsigned int *subtree_member_list_indices);

#endif
