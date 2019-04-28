#include<cstdlib>
#include<cstring>
#include<cmath>
#include<cctype>

#include<iostream>
#include<vector>
#include<map>
#include <chrono>

#include "single_tree_helper.h"

using namespace std;
using namespace std::chrono;

void run_benchmark(unsigned int tree_depth) {
    unsigned long long tree_size = 1 << tree_depth;
    //unsigned int log_num_labels = 3 * tree_depth / 4;
    unsigned int log_num_labels = tree_depth - 4;
    unsigned long long num_labels = 1 << log_num_labels;
    unsigned int subtree_depth = tree_depth - log_num_labels;
    unsigned int threshold = tree_depth;

    unsigned long long depths[num_labels];
    unsigned long long labels[num_labels];

    for (unsigned long long i = 0; i < num_labels; ++i) {
        unsigned long long l = 1 + (unsigned long long)(rand() % (tree_size - 1));
        depths[i] = floor(log2(l));
        labels[i] = l;
        //printf("Label: %d\t Depth: %d\n", labels[i], depths[i]);
    }

    printf("Start recording time.\n\n");
    auto start = high_resolution_clock::now();

    pair< MapLabelToLabelVector, vector<unsigned long long> > X = create_subtree_partitions(labels, depths, num_labels, tree_depth, subtree_depth);
    MapLabelToLabelVector subtree_label_to_label_index_list = X.first;
    vector<unsigned long long> subtree_labels = X.second;

    vector< pair<unsigned long long, unsigned long long> > edges = get_all_edges(
            labels, depths, threshold, tree_depth, subtree_depth, subtree_labels, subtree_label_to_label_index_list);

    auto stop = high_resolution_clock::now();
    int64_t elapsed_time = duration_cast<nanoseconds>(stop - start).count();
    printf("Elapsed Time: %lf miliseconds\n", elapsed_time / 1e6);
    printf("Found %ld edges amongst %lld vertices.\n", edges.size(), num_labels);
    printf("%lf nanoseconds per edge\n", 1.0 * elapsed_time / edges.size());

}


int main(int argc, char *argv[]) {

    if (argc < 2) {
        printf("Not enough arguments!\n");
        return 1;
    }
    else if (argc > 2) {
        printf("Too many arguments!\n");
        return 1;
    }
    else {
        for (unsigned long long i = 0; i < strlen(argv[1]); ++i) {
            if (! isdigit(argv[1][i])) {
                printf("Invalid Argument! Expected \"unsigned int\". Received: %s\n", argv[1]);
                return 1;
            }
        }
    }

    unsigned int tree_depth = (unsigned int)(stoi(argv[1]));

    run_benchmark(tree_depth);

    return 0;
}
