#include "product_trees_helper.h"

#define ROOT_SUBTREE_LABEL 1

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
        depth = floor(log2(label));
        extra_depth = depth % subtree_depth;
        subtree_label_to_label_index_list[label >> extra_depth].push_back(label);
    }
    //for (auto& t : subtree_label_to_label_index_list) {
        //cout << t.first << "\t";
        //for (const auto i: t.second)
            //cout << i << " ";
        //cout << endl;
    //}

    printf("Found %ld subtrees.\n", subtree_labels.size());
    return pair<MapLabelToLabelVector, vector<unsigned long long> >(subtree_label_to_label_index_list, subtree_labels);
}

void add_descendant_subtree_labels_to_vector(
        unsigned long long source_subtree_label, unsigned int tree_depth, unsigned int subtree_depth,
        unsigned int remaining_distance, vector<unsigned long long> &subtree_labels_vector) {
    unsigned int current_depth = floor(log2(source_subtree_label));
    unsigned long long partial_label_max = 1;
    unsigned long long current_base_label = source_subtree_label;
    while (current_depth < tree_depth && remaining_distance >= 0) {
        for (unsigned long long partial_label = 0; partial_label < partial_label_max; ++partial_label) {
            unsigned long long descendant_subtree_label = current_base_label + partial_label;
            subtree_labels_vector.push_back(descendant_subtree_label);
        }
        current_base_label = current_base_label << subtree_depth;
        current_depth += subtree_depth;
        partial_label_max = partial_label_max << subtree_depth;
        remaining_distance -= subtree_depth;
    }
}

void get_nearby_subtrees(
        unsigned long long source_subtree_label, unsigned int tree_depth, unsigned int subtree_depth,
        unsigned int remaining_distance, vector<unsigned long long> &visited_subtree_labels) {
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

    unsigned long long current_root_subtree_label = source_subtree_label;
    unsigned long long old_root_subtree_label;
    unsigned long long neighboring_subtree_label;
    unsigned long long subtree_valence = 1 << subtree_depth;

    add_descendant_subtree_labels_to_vector(source_subtree_label, tree_depth, subtree_depth, remaining_distance, visited_subtree_labels);
    while (remaining_distance >= subtree_depth) {
        /* TODO: generate all subtrees below current one <30-04-19, shankha> */
        if (current_root_subtree_label != ROOT_SUBTREE_LABEL) {
            old_root_subtree_label = current_root_subtree_label;
            current_root_subtree_label = current_root_subtree_label >> subtree_depth;
            remaining_distance -= 2 * subtree_depth;
            for (unsigned long long partial_label = 0; partial_label < subtree_valence; ++partial_label) {
                neighboring_subtree_label = (current_root_subtree_label << subtree_depth) + partial_label;
                add_descendant_subtree_labels_to_vector(
                        neighboring_subtree_label, tree_depth, subtree_depth, remaining_distance, visited_subtree_labels);
            }
            remaining_distance += subtree_depth;
        }
        else
            break;
    }

}


vector< pair<unsigned long long, unsigned long long> > get_all_edges(
        unsigned long long *labels, unsigned int threshold,
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
        first_subtree_member_labels = subtree_label_to_label_index_list[first_subtree_label];
        visited_subtree_labels.clear();
        unsigned int remaining_distance = threshold;
        get_nearby_subtrees(first_subtree_label, tree_depth, subtree_depth, remaining_distance, visited_subtree_labels);
        for (unsigned long long second_subtree_label : visited_subtree_labels) {
            if (second_subtree_label >= first_subtree_label) {

        //for (unsigned long long j = i + 1; j < num_subtrees; ++j) {
            //unsigned long long second_subtree_label = subtree_labels[j];
            //unsigned int dist = label_distance(first_subtree_label, second_subtree_label);

            //if (dist <= threshold) {
                second_subtree_member_labels = subtree_label_to_label_index_list[second_subtree_label];
                degree += first_subtree_member_labels.size() * second_subtree_member_labels.size();

                for (unsigned int k = 0; k < first_subtree_member_labels.size(); ++k) {
                    unsigned long long first_endpoint_label = first_subtree_member_labels[k];
                    for (unsigned int l = 0; l < second_subtree_member_labels.size(); ++l) {
                        unsigned long long second_endpoint_label = second_subtree_member_labels[l];
                        if (label_distance(first_endpoint_label, second_endpoint_label) <= threshold)
                            edges.push_back(pair<unsigned long long, unsigned long long>(first_endpoint_label, second_endpoint_label));
                    }
                }
            }
        }

        for (unsigned int k = 0; k < first_subtree_member_labels.size(); ++k) {
            unsigned long long first_endpoint_label = first_subtree_member_labels[k];
            for (unsigned int l = k + 1; l < first_subtree_member_labels.size(); ++l) {
                unsigned long long second_endpoint_label = first_subtree_member_labels[l];
                if (label_distance(first_endpoint_label, second_endpoint_label) <= threshold)
                    edges.push_back(pair<unsigned long long, unsigned long long>(first_endpoint_label, second_endpoint_label));
            }
        }
    }

    return edges;
}
