#include<cstdlib>
#include<cstring>
#include<cmath>
#include<cctype>

#include<iostream>
#include<vector>
#include<map>
#include <chrono>

#include "product_trees_helper.h"

using namespace std;
using namespace std::chrono;

void run_benchmark(unsigned int tree_depth, unsigned int subtree_depth, unsigned int log_num_labels) {
    unsigned int num_trees = 1;

    unsigned long long tree_size = 1 << tree_depth;
    unsigned long long num_labels = 1 << log_num_labels;
    unsigned int threshold = num_trees * tree_depth;

    unsigned long long labels[num_labels];

    for (unsigned long long i = 0; i < num_labels; ++i) {
        unsigned long long l = 1 + (unsigned long long)(rand() % (tree_size - 1));
        labels[i] = l;
    }

    for (unsigned int i = 0; i < num_labels; ++i) {
        //printf("%lld\n", labels[i]);
    }

    printf("Start recording time.\n\n");
    auto start = high_resolution_clock::now();

    vector<unsigned long long> subtree_labels;
    unordered_map<unsigned long long, unsigned int> map_subtree_label_to_index;
    unsigned long long *subtree_label_to_label_index_list;
    unsigned int *subtree_member_list_indices;

    tie(subtree_labels, map_subtree_label_to_index) = get_subtree_labels(labels, num_labels, tree_depth, subtree_depth);
    tie(subtree_label_to_label_index_list, subtree_member_list_indices) = create_subtree_partitions(
            labels, num_labels, tree_depth, subtree_depth, subtree_labels, map_subtree_label_to_index);


    vector< pair<unsigned long long, unsigned long long> > edges = get_all_edges(
            labels, threshold, tree_depth, subtree_depth,
            subtree_labels, map_subtree_label_to_index,
            subtree_label_to_label_index_list, subtree_member_list_indices);

    auto stop = high_resolution_clock::now();
    int64_t elapsed_time = duration_cast<nanoseconds>(stop - start).count();
    printf("Elapsed Time: %lf miliseconds\n", elapsed_time / 1e6);
    printf("Found %ld edges amongst %lld vertices.\n", edges.size(), num_labels);
    printf("%lf nanoseconds per edge\n", 1.0 * elapsed_time / edges.size());

    for (pair<unsigned long long, unsigned long long> edge: edges) {
        //printf("%lld\t%lld\n", edge.first, edge.second);
    }

}


int main(int argc, char *argv[]) {

    if (argc < 3) {
        printf("Not enough arguments!\n");
        return 1;
    }
    else if (argc > 4) {
        printf("Too many arguments!\n");
        return 1;
    }
    else {
        for (unsigned long long argv_index = 1; argv_index < 4; ++argv_index) {
            for (unsigned long long i = 0; i < strlen(argv[argv_index]); ++i) {
                if (! isdigit(argv[argv_index][i])) {
                    printf("Invalid Argument! Expected \"unsigned int\". Received: \"%s\"\n", argv[1]);
                    return 1;
                }
            }
        }
    }

    unsigned int log_num_labels = (unsigned int)(stoi(argv[1]));
    unsigned int tree_depth = (unsigned int)(stoi(argv[2]));
    unsigned int subtree_depth = (unsigned int)(stoi(argv[3]));

    run_benchmark(tree_depth, subtree_depth, log_num_labels);

    return 0;
}
