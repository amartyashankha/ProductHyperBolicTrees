#include "single_tree_helper.h"

#define ROOT_SUBTREE_LABEL 1

inline unsigned int label_distance(LONG_UINT first_label, LONG_UINT second_label) {
    unsigned int dist = 0;
    LONG_UINT l1, l2;
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


tuple<vector<LONG_UINT>, UNORDERED_MAP> get_subtree_labels(
        LONG_UINT *labels, LONG_UINT num_labels, unsigned int tree_depth, unsigned int subtree_depth) {
    UNORDERED_MAP map_subtree_label_to_index;
    vector<LONG_UINT> subtree_labels;
    LONG_UINT subtree_label;

    unsigned int current_depth = 0;
    while (current_depth < tree_depth) {
        LONG_UINT current_base_label = (1 << current_depth);
        for (LONG_UINT partial_label = 0; partial_label < current_base_label; ++partial_label) {
            subtree_label = partial_label + current_base_label;
            subtree_labels.push_back(subtree_label);
            map_subtree_label_to_index[subtree_label] = subtree_labels.size() - 1;
        }
        current_depth += subtree_depth;
    }
    printf("Found %ld subtrees.\n", subtree_labels.size());

    return make_tuple(subtree_labels, map_subtree_label_to_index);
}


tuple< LONG_UINT *, unsigned int * > create_subtree_partitions(
        LONG_UINT *labels, LONG_UINT num_labels,
        unsigned int tree_depth, unsigned int subtree_depth,
        vector<LONG_UINT> subtree_labels, UNORDERED_MAP map_subtree_label_to_index) {
    unsigned int extra_depth, subtree_index;
    LONG_UINT subtree_label, depth;

    LONG_UINT *subtree_label_to_label_index_list = (LONG_UINT *)malloc(num_labels * sizeof(LONG_UINT));
    unsigned int num_subtrees = subtree_labels.size();
    unsigned int subtree_member_list_sizes[num_subtrees];
    unsigned int subtree_member_list_pointers[num_subtrees + 1];
    unsigned int *subtree_member_list_indices = (unsigned int *)malloc((num_subtrees + 1) * sizeof(unsigned int));

    for (unsigned int i = 0; i < num_subtrees; ++i)
        subtree_member_list_sizes[i] = 0;
    for (LONG_UINT i = 0; i < num_labels; ++i) {
        subtree_label = labels[i];
        depth = floor(log2(subtree_label));
        extra_depth = depth % subtree_depth;
        subtree_label = subtree_label >> extra_depth;
        ++subtree_member_list_sizes[map_subtree_label_to_index[subtree_label]];
    }

    // Compute prefix sums
    subtree_member_list_pointers[0] = 0;
    subtree_member_list_indices[0] = 0;
    for (unsigned int i = 1; i < num_subtrees + 1; ++i) {
        subtree_member_list_pointers[i] = subtree_member_list_pointers[i - 1] + subtree_member_list_sizes[i - 1];
        subtree_member_list_indices[i] = subtree_member_list_pointers[i];
    }

    for (LONG_UINT i = 0; i < num_labels; ++i) {
        LONG_UINT label = labels[i];
        depth = floor(log2(label));
        extra_depth = depth % subtree_depth;
        subtree_label = label >> extra_depth;
        subtree_index = map_subtree_label_to_index[subtree_label];
        subtree_label_to_label_index_list[subtree_member_list_pointers[subtree_index]++] = label;
    }

    //for (unsigned int i = 0; i < num_subtrees; ++i)
        //cout << subtree_labels[i] << "\t" << subtree_member_list_sizes[i] << " "<< subtree_member_list_indices[i] << endl;
    //cout << endl;
    //for (unsigned int i = 0; i < num_subtrees; ++i) {
        //cout << subtree_labels[i] << "\t";
        //for (unsigned int j = subtree_member_list_indices[i]; j < subtree_member_list_indices[i + 1]; ++j)
            //cout << subtree_label_to_label_index_list[j] << " ";
        //cout << endl;
    //}
    //cout << endl;

    return tuple< LONG_UINT *, unsigned int * >(subtree_label_to_label_index_list, subtree_member_list_indices);
}

void add_descendant_subtree_labels_to_vector(
        LONG_UINT source_subtree_label, unsigned int tree_depth, unsigned int subtree_depth,
        int remaining_distance, vector<LONG_UINT> &subtree_labels_vector) {
    unsigned int current_depth = floor(log2(source_subtree_label));
    LONG_UINT partial_label_max = 1;
    LONG_UINT current_base_label = source_subtree_label;

    while (current_depth < tree_depth && remaining_distance >= 0) {
        for (LONG_UINT partial_label = 0; partial_label < partial_label_max; ++partial_label) {
            LONG_UINT descendant_subtree_label = current_base_label + partial_label;
            subtree_labels_vector.push_back(descendant_subtree_label);
        }
        current_base_label = current_base_label << subtree_depth;
        current_depth += subtree_depth;
        partial_label_max = partial_label_max << subtree_depth;
        remaining_distance -= subtree_depth;
    }
}

void get_nearby_subtrees(
        LONG_UINT source_subtree_label, unsigned int tree_depth, unsigned int subtree_depth,
        int remaining_distance, vector<LONG_UINT> &visited_subtree_labels) {
    LONG_UINT current_root_subtree_label = source_subtree_label;
    LONG_UINT old_root_subtree_label;
    LONG_UINT neighboring_subtree_label;
    LONG_UINT subtree_valence = 1 << subtree_depth;
    unsigned int current_remaining_distance = remaining_distance;

    add_descendant_subtree_labels_to_vector(source_subtree_label, tree_depth, subtree_depth, remaining_distance, visited_subtree_labels);

    while (current_remaining_distance >= subtree_depth) {
        if (current_root_subtree_label != ROOT_SUBTREE_LABEL) {
            old_root_subtree_label = current_root_subtree_label;
            current_root_subtree_label = current_root_subtree_label >> subtree_depth;
            visited_subtree_labels.push_back(current_root_subtree_label);
            for (LONG_UINT partial_label = 0; partial_label < subtree_valence; ++partial_label) {
                neighboring_subtree_label = (current_root_subtree_label << subtree_depth) + partial_label;
                unsigned int recursive_remaining_distance = remaining_distance - label_distance(source_subtree_label, neighboring_subtree_label);
                if (neighboring_subtree_label != old_root_subtree_label)
                    add_descendant_subtree_labels_to_vector(
                            neighboring_subtree_label, tree_depth, subtree_depth, recursive_remaining_distance, visited_subtree_labels);
            }
            current_remaining_distance -= subtree_depth;
        }
        else
            break;
    }
    //cout << source_subtree_label << "\t";
    //for (auto neighboring_subtree_label : visited_subtree_labels)
        //cout << " " << neighboring_subtree_label;
    //cout << endl;

}


vector< pair<LONG_UINT, LONG_UINT> > get_all_edges(
        LONG_UINT *labels, LONG_UINT num_labels, unsigned int threshold, unsigned int tree_depth, unsigned int subtree_depth,
        vector<LONG_UINT> subtree_labels, UNORDERED_MAP map_subtree_label_to_index,
        LONG_UINT *subtree_label_to_label_index_list, unsigned int *subtree_member_list_indices) {
    vector< pair<LONG_UINT, LONG_UINT> > edges;
    edges.reserve(100 * num_labels);

    LONG_UINT num_subtrees = subtree_labels.size();
    LONG_UINT *first_subtree_member_labels, *second_subtree_member_labels;
    unsigned int num_first_subtree_member_labels, num_second_subtree_member_labels;
    vector<LONG_UINT> visited_subtree_labels;
    vector<unsigned int> subtree_neighbors_accumulator;

    for (LONG_UINT i = 0; i < num_subtrees; ++i) {
        subtree_neighbors_accumulator.push_back(0);
        LONG_UINT first_subtree_label = subtree_labels[i];
        first_subtree_member_labels = &subtree_label_to_label_index_list[subtree_member_list_indices[i]];
        num_first_subtree_member_labels = subtree_member_list_indices[i + 1] - subtree_member_list_indices[i];

        visited_subtree_labels.clear();
        int remaining_distance = threshold + 0 * subtree_depth;
        get_nearby_subtrees(first_subtree_label, tree_depth, subtree_depth, remaining_distance, visited_subtree_labels);
        for (LONG_UINT second_subtree_label : visited_subtree_labels) {
            ++subtree_neighbors_accumulator[i];
            if (second_subtree_label > first_subtree_label) {
                unsigned int second_subtree_index = map_subtree_label_to_index[second_subtree_label];
                second_subtree_member_labels = &subtree_label_to_label_index_list[subtree_member_list_indices[second_subtree_index]];
                num_second_subtree_member_labels = subtree_member_list_indices[second_subtree_index + 1];
                num_second_subtree_member_labels -= subtree_member_list_indices[second_subtree_index];

                for (unsigned int k = 0; k < num_first_subtree_member_labels; ++k) {
                    LONG_UINT first_endpoint_label = first_subtree_member_labels[k];
                    for (unsigned int l = 0; l < num_second_subtree_member_labels; ++l) {
                        LONG_UINT second_endpoint_label = second_subtree_member_labels[l];
                        if (label_distance(first_endpoint_label, second_endpoint_label) <= threshold)
                            edges.push_back(pair<LONG_UINT, LONG_UINT>(first_endpoint_label, second_endpoint_label));
                    }
                }
            }
        }

        for (unsigned int k = 0; k < num_first_subtree_member_labels; ++k) {
            LONG_UINT first_endpoint_label = first_subtree_member_labels[k];
            for (unsigned int l = k + 1; l < num_first_subtree_member_labels; ++l) {
                LONG_UINT second_endpoint_label = first_subtree_member_labels[l];
                if (label_distance(first_endpoint_label, second_endpoint_label) <= threshold)
                    edges.push_back(pair<LONG_UINT, LONG_UINT>(first_endpoint_label, second_endpoint_label));
            }
        }
    }
    cout << "Average subtree neighborhood size: "
         << accumulate(subtree_neighbors_accumulator.begin(), subtree_neighbors_accumulator.end(), 0) * 1.0 / num_subtrees
         << endl << endl;

    return edges;
}
