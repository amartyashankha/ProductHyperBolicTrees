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

void run_benchmark(unsigned int tree_depth, unsigned int subtree_depth) {
    unsigned int num_trees = 1;

    unsigned long long tree_size = 1 << tree_depth;
    unsigned int log_num_labels = 10;
    unsigned long long num_labels = 1 << log_num_labels;
    unsigned int threshold = num_trees * tree_depth;

    unsigned long long depths[num_labels];
    unsigned long long labels[num_labels];

    for (unsigned long long i = 0; i < num_labels; ++i) {
        unsigned long long l = 1 + (unsigned long long)(rand() % (tree_size - 1));
        depths[i] = floor(log2(l));
        labels[i] = l;
    }

    for (int i = 0; i < num_labels; ++i) {
        //printf("%lld\n", labels[i]);
    }

    printf("Start recording time.\n\n");
    auto start = high_resolution_clock::now();

    pair< MapLabelToLabelVector, vector<unsigned long long> > X = create_subtree_partitions(labels, depths, num_labels, tree_depth, subtree_depth);
    MapLabelToLabelVector subtree_label_to_label_index_list = X.first;
    vector<unsigned long long> subtree_labels = X.second;

    vector< pair<unsigned long long, unsigned long long> > edges = get_all_edges(
            labels, threshold, tree_depth, subtree_depth, subtree_labels, subtree_label_to_label_index_list);

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
    else if (argc > 3) {
        printf("Too many arguments!\n");
        return 1;
    }
    else {
        for (unsigned long long argv_index = 1; argv_index < 3; ++argv_index) {
            for (unsigned long long i = 0; i < strlen(argv[argv_index]); ++i) {
                if (! isdigit(argv[1][i])) {
                    printf("Invalid Argument! Expected \"unsigned int\". Received: %s\n", argv[1]);
                    return 1;
                }
            }
        }
    }

    unsigned int tree_depth = (unsigned int)(stoi(argv[1]));
    unsigned int subtree_depth = (unsigned int)(stoi(argv[2]));

    run_benchmark(tree_depth, subtree_depth);

    return 0;
}
