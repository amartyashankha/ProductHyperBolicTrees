#ifndef PRODUCT_TREES_HELPER_H
#define PRODUCT_TREES_HELPER_H

#include<iostream>
#include<vector>
#include<numeric>
#include<unordered_map>
#include<unordered_set>

#include<cassert>
#include<cmath>

#include <sparsehash/dense_hash_set>

using namespace std;

struct pair_hash {
    template <class T1, class T2>
    std::size_t operator () (const std::pair<T1,T2> &p) const {
        auto h1 = std::hash<T1>{}(p.first);
        auto h2 = std::hash<T2>{}(p.second);

        // Mainly for demonstration purposes, i.e. works but is overly simple
        // In the real world, use sth. like boost.hash_combine
        return h1 ^ h2;
    }
};

typedef unsigned int LONG_UINT;
typedef pair<LONG_UINT, LONG_UINT> LABEL;
typedef unordered_map<LABEL, unsigned int, pair_hash> UNORDERED_MAP;
typedef unordered_set<LABEL, pair_hash> UNORDERED_SET;

tuple<vector<LABEL>, UNORDERED_MAP> get_subtree_labels(
        LABEL *labels, LONG_UINT num_labels, unsigned int tree_depth, unsigned int subtree_depth);

tuple< LABEL *, unsigned int * > create_subtree_partitions(
        LABEL *labels, LONG_UINT num_labels,
        unsigned int tree_depth, unsigned int subtree_depth,
        vector<LABEL> subtree_labels, UNORDERED_MAP map_subtree_label_to_index);

vector< pair<LABEL, LABEL> > get_all_edges(
        LABEL *labels, LONG_UINT num_labels, unsigned int threshold, unsigned int tree_depth, unsigned int subtree_depth,
        vector<LABEL> subtree_labels, UNORDERED_MAP map_subtree_label_to_index,
        LABEL *subtree_label_to_label_index_list, unsigned int *subtree_member_list_indices);

#endif
