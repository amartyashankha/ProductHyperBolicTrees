#include "product_trees_helper.h"


unsigned int label_distance(unsigned long long first_label, unsigned long long second_label) {
    unsigned int dist = 0;
    unsigned long long l1, l2;
    unsigned int height_delta;
    unsigned int first_label_depth = 32 - __builtin_clz(first_label);
    unsigned int second_label_depth = 32 - __builtin_clz(second_label);

    if (second_label_depth > first_label_depth) {
        height_delta = second_label_depth - first_label_depth;
        l1 = first_label;
        l2 = second_label;
    }
    else {
        height_delta = first_label_depth - second_label_depth;
        l1 = second_label;
        l2 = first_label;
    }
    l2 = l2 >> height_delta;
    dist = 2 * (31 ^ __builtin_clz(((l2 ^ l1) << 1) + 1)) + height_delta;

    return dist;
}


pair<MapLabelToLabelVector, vector<unsigned long long> > create_subtree_partitions(
        unsigned long long *labels, unsigned long long *depths, unsigned long long num_labels,
        unsigned int tree_depth, unsigned int subtree_depth) {
    unsigned int extra_depth;
    unsigned long long label, depth;

    MapLabelToLabelVector subtree_label_to_label_index_list;
    vector<unsigned long long> empty_vector;
    vector<unsigned long long> subtree_labels;

    unsigned int current_depth = 0;
    while (current_depth < tree_depth) {
        unsigned long long current_base_label = (1 << current_depth);
        for (unsigned long long partial_label = 0; partial_label < current_base_label; ++partial_label) {
            label = partial_label + current_base_label;
            subtree_label_to_label_index_list[label] = empty_vector;
            subtree_labels.push_back(label);
        }
        current_depth += subtree_depth;
    }

    for (unsigned long long i = 0; i < num_labels; ++i) {
        label = labels[i];
        depth = depths[i];
        extra_depth = depth % subtree_depth;
        subtree_label_to_label_index_list[label >> extra_depth].push_back(label);
    }

    printf("Found %ld subtrees.\n", subtree_labels.size());
    return pair<MapLabelToLabelVector, vector<unsigned long long> >(subtree_label_to_label_index_list, subtree_labels);
}

void get_nearby_subtrees(
        unsigned long long subtree_label, unsigned int tree_depth, unsigned int subtree_depth,
        unsigned int remaining_distance, unordered_set<unsigned long long> &visited_subtree_labels) {
    unsigned int subtree_label_depth = floor(log2(subtree_label));

    /* TODO: Remove this for optimization later. <16-04-19, shankha> */
    //if ((subtree_label_depth % subtree_depth) != 0) {
        //printf("Invalid subtree label: %lld. Label depth: %d must be a multiple of the subtree depth: %d.\n",
                //subtree_label, subtree_label_depth, subtree_depth);
        //assert(0);
    //}
    //if (subtree_label_depth > tree_depth - subtree_depth) {
        //printf("Invalid subtree label: %lld. Label depth: %d must be a less than: %d.\n",
                //subtree_label, subtree_label_depth, tree_depth - subtree_depth);
        //assert(0);
    //}
    //if ((remaining_distance % subtree_depth) != 0) {
        //printf("Invalid remaining distance: %d. Must be a multiple of the subtree depth: %d.\n", remaining_distance, subtree_depth);
        //assert(0);
    //}

    if (visited_subtree_labels.find(subtree_label) != visited_subtree_labels.end())
        return;

    visited_subtree_labels.insert(subtree_label);
    if (remaining_distance < subtree_depth)
        return;
    else {

        /* TODO: This is not quite correct? <16-04-19, shankha> */
        remaining_distance -= subtree_depth;
        if (subtree_label_depth + subtree_depth <= tree_depth - subtree_depth) {
            unsigned long long child_subtree_label_base = subtree_label << subtree_depth;
            unsigned long long max_label_extension = 1 << subtree_depth;
            for (unsigned long long label_extension = 0; label_extension < max_label_extension; ++label_extension) {
                unsigned long long child_subtree_label = child_subtree_label_base + label_extension;
                get_nearby_subtrees(child_subtree_label, tree_depth, subtree_depth, remaining_distance, visited_subtree_labels);
            }
        }

        if (subtree_label_depth >= subtree_depth) {
            unsigned long long parent_subtree_label = subtree_label >> subtree_depth;
            get_nearby_subtrees(parent_subtree_label, tree_depth, subtree_depth, remaining_distance, visited_subtree_labels);
        }
    }
}


vector< pair<unsigned long long, unsigned long long> > get_all_edges(
        unsigned long long *labels, unsigned long long *depths, unsigned int threshold,
        unsigned int tree_depth, unsigned int subtree_depth,
        vector<unsigned long long> subtree_labels,
        MapLabelToLabelVector subtree_label_to_label_index_list) {
    unsigned int degree = 0;
    vector< pair<unsigned long long, unsigned long long> > edges;

    //pair< MapLabelToLabelVector, vector<unsigned long long> > X = create_subtree_partitions(labels, depths, num_labels, tree_depth, subtree_depth);
    //MapLabelToLabelVector subtree_label_to_label_index_list = X.first;
    //vector<unsigned long long> subtree_labels = X.second;
    unsigned long long num_subtrees = subtree_labels.size();
    vector<unsigned long long> first_subtree_member_labels, second_subtree_member_labels;
    unordered_set<unsigned long long> visited_subtree_labels;

    for (unsigned long long i = 0; i < num_subtrees; ++i) {
        unsigned long long first_subtree_label = subtree_labels[i];
        visited_subtree_labels.clear();
        get_nearby_subtrees(first_subtree_label, tree_depth, subtree_depth, threshold, visited_subtree_labels);
        for (unsigned long long second_subtree_label : visited_subtree_labels) {
            if (second_subtree_label >= first_subtree_label) {

        //for (unsigned long long j = i; j < num_subtrees; ++j) {
            //unsigned long long second_subtree_label = subtree_labels[j];
            //unsigned int dist = label_distance(first_subtree_label, second_subtree_label);

            //if (dist <= threshold) {
                first_subtree_member_labels = subtree_label_to_label_index_list[first_subtree_label];
                second_subtree_member_labels = subtree_label_to_label_index_list[second_subtree_label];
                degree += first_subtree_member_labels.size() * second_subtree_member_labels.size();

                for (unsigned int k = 0; k < first_subtree_member_labels.size(); ++k) {
                    unsigned long long first_endpoint_label = first_subtree_member_labels[k];
                    for (unsigned int l = 0; l < first_subtree_member_labels.size(); ++l) {
                        unsigned long long second_endpoint_label = second_subtree_member_labels[k];
                        if(label_distance(first_endpoint_label, second_endpoint_label) <= threshold)
                            edges.push_back(pair<unsigned long long, unsigned long long>(first_endpoint_label, second_endpoint_label));
                    }
                }
                // print(dist, first_subtree_member_labels.size(), second_subtree_member_labels.size())
            }
        }
    }

    return edges;
}
