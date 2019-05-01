#ifndef PRODUCT_TREES_HELPER_H
#define PRODUCT_TREES_HELPER_H

#include<iostream>
#include<vector>
#include<numeric>
#include<unordered_map>

#include<cassert>
#include<cmath>

using namespace std;

typedef unordered_map< unsigned long long, vector<unsigned long long> > MapLabelToLabelVector;

tuple< vector<unsigned long long>, unordered_map<unsigned long long, unsigned int> > get_subtree_labels(
        unsigned long long *labels, unsigned long long num_labels, unsigned int tree_depth, unsigned int subtree_depth);

tuple< unsigned long long *, unsigned int * > create_subtree_partitions(
        unsigned long long *labels, unsigned long long num_labels,
        unsigned int tree_depth, unsigned int subtree_depth,
        vector<unsigned long long> subtree_labels, unordered_map<unsigned long long, unsigned int> map_subtree_label_to_index);

vector< pair<unsigned long long, unsigned long long> > get_all_edges(
        unsigned long long *labels, unsigned int threshold, unsigned int tree_depth, unsigned int subtree_depth,
        vector<unsigned long long> subtree_labels, unordered_map<unsigned long long, unsigned int> map_subtree_label_to_index,
        unsigned long long *subtree_label_to_label_index_list, unsigned int *subtree_member_list_indices);

#endif
