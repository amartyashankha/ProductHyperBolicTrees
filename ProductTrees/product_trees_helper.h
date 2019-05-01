#ifndef PRODUCT_TREES_HELPER_H
#define PRODUCT_TREES_HELPER_H

#include<iostream>
#include<vector>
#include<numeric>
#include<unordered_map>

#include<cassert>
#include<cmath>

#include "absl/container/flat_hash_map.h"

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

LABEL operator %(const LABEL &p, const unsigned int &a) {
    return make_pair(p.first % a, p.second % a);
}

LABEL operator >>(const LABEL &p, const LABEL &a) {
    return make_pair(p.first >> a.first, p.second >> a.second);
}

tuple<vector<LABEL>, UNORDERED_MAP> get_subtree_labels(
        LABEL *labels, LONG_UINT num_labels, unsigned int tree_depth, unsigned int subtree_depth);

tuple< LABEL *, unsigned int * > create_subtree_partitions(
        LABEL *labels, LONG_UINT num_labels,
        unsigned int tree_depth, unsigned int subtree_depth,
        vector<LABEL> subtree_labels, UNORDERED_MAP map_subtree_label_to_index);

vector< pair<LONG_UINT, LONG_UINT> > get_all_edges(
        LONG_UINT *labels, unsigned int threshold, unsigned int tree_depth, unsigned int subtree_depth,
        vector<LONG_UINT> subtree_labels, UNORDERED_MAP map_subtree_label_to_index,
        LONG_UINT *subtree_label_to_label_index_list, unsigned int *subtree_member_list_indices);

#endif
