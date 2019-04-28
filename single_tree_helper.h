#ifndef SINGLE_TREE_HELPER_H
#define SINGLE_TREE_HELPER_H

#include<iostream>
#include<vector>
#include<unordered_map>
#include<unordered_set>

#include<cassert>
#include<cmath>

using namespace std;

typedef unordered_map< unsigned long long, vector<unsigned long long> > MapLabelToLabelVector;

pair<MapLabelToLabelVector, vector<unsigned long long> > create_subtree_partitions(
        unsigned long long *labels, unsigned long long *depths, unsigned long long num_labels,
        unsigned int tree_depth, unsigned int subtree_depth);

vector< pair<unsigned long long, unsigned long long> > get_all_edges(
        unsigned long long *labels, unsigned long long *depths, unsigned int threshold,
        unsigned int tree_depth, unsigned int subtree_depth,
        vector<unsigned long long> subtree_labels,
        MapLabelToLabelVector subtree_label_to_label_index_list);

#endif
