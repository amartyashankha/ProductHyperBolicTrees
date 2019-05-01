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


tuple< vector<unsigned long long>, unordered_map<unsigned long long, unsigned int> > get_subtree_labels(
        unsigned long long *labels, unsigned long long num_labels, unsigned int tree_depth, unsigned int subtree_depth) {
    unordered_map<unsigned long long, unsigned int> map_subtree_label_to_index;
    vector<unsigned long long> subtree_labels;
    unsigned long long subtree_label;

    unsigned int current_depth = 0;
    while (current_depth < tree_depth) {
        unsigned long long current_base_label = (1 << current_depth);
        for (unsigned long long partial_label = 0; partial_label < current_base_label; ++partial_label) {
            subtree_label = partial_label + current_base_label;
            subtree_labels.push_back(subtree_label);
            map_subtree_label_to_index[subtree_label] = subtree_labels.size() - 1;
        }
        current_depth += subtree_depth;
    }
    printf("Found %ld subtrees.\n", subtree_labels.size());

    return make_tuple(subtree_labels, map_subtree_label_to_index);
}


tuple< unsigned long long *, unsigned int * > create_subtree_partitions(
        unsigned long long *labels, unsigned long long num_labels,
        unsigned int tree_depth, unsigned int subtree_depth,
        vector<unsigned long long> subtree_labels, unordered_map<unsigned long long, unsigned int> map_subtree_label_to_index) {
    unsigned int extra_depth, subtree_index;
    unsigned long long subtree_label, depth;

    unsigned long long *subtree_label_to_label_index_list = (unsigned long long *)malloc(num_labels * sizeof(unsigned long long));
    unsigned int num_subtrees = subtree_labels.size();
    unsigned int subtree_member_list_sizes[num_subtrees];
    unsigned int subtree_member_list_pointers[num_subtrees + 1];
    unsigned int *subtree_member_list_indices = (unsigned int *)malloc((num_subtrees + 1) * sizeof(unsigned int));

    for (unsigned int i = 0; i < num_subtrees; ++i)
        subtree_member_list_sizes[i] = 0;
    for (unsigned long long i = 0; i < num_labels; ++i) {
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

    for (unsigned long long i = 0; i < num_labels; ++i) {
        unsigned long long label = labels[i];
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

    return tuple< unsigned long long *, unsigned int * >(subtree_label_to_label_index_list, subtree_member_list_indices);
}

void add_descendant_subtree_labels_to_vector(
        unsigned long long source_subtree_label, unsigned int tree_depth, unsigned int subtree_depth,
        unsigned int max_depth, unsigned long long original_subtree_label,
        int remaining_distance, vector<unsigned long long> &subtree_labels_vector) {
    unsigned int source_depth = floor(log2(source_subtree_label));
    unsigned int current_depth = source_depth;
    unsigned long long partial_label_max = 1;
    unsigned long long current_base_label = source_subtree_label;
    //cout << "Down from: " << source_subtree_label << " Remaining distance: " << remaining_distance << endl;

    while (current_depth <= max_depth && remaining_distance >= 0) {
        if (current_depth == max_depth) {
            if (current_base_label > original_subtree_label) {
                //unsigned int extra_depth = current_depth - source_depth;
                //partial_label_max = ((1 << extra_depth) - 1) & original_subtree_label;
                partial_label_max = 0;
                //cout << endl << "\tPartial labels: " << current_base_label << " " << partial_label_max << endl;
            }
        }
        for (unsigned long long partial_label = 0; partial_label < partial_label_max; ++partial_label) {
            unsigned long long descendant_subtree_label = current_base_label + partial_label;
            subtree_labels_vector.push_back(descendant_subtree_label);
            //cout << descendant_subtree_label << " ";
        }
        current_base_label = current_base_label << subtree_depth;
        current_depth += subtree_depth;
        partial_label_max = partial_label_max << subtree_depth;
        remaining_distance -= subtree_depth;
    }
    //cout << endl;
}

void get_nearby_subtrees(
        unsigned long long source_subtree_label, unsigned int tree_depth, unsigned int subtree_depth,
        int remaining_distance, vector<unsigned long long> &visited_subtree_labels) {
    unsigned long long current_root_subtree_label = source_subtree_label;
    unsigned long long old_root_subtree_label;
    unsigned long long neighboring_subtree_label;
    unsigned long long subtree_valence = 1 << subtree_depth;
    unsigned int current_remaining_distance = remaining_distance;
    unsigned int max_depth = floor(log2(source_subtree_label));
    //cout << "Processing subtree " << source_subtree_label << endl;

    add_descendant_subtree_labels_to_vector(
            source_subtree_label, tree_depth, subtree_depth, max_depth, source_subtree_label, remaining_distance, visited_subtree_labels);

    while (current_remaining_distance >= subtree_depth) {
        if (current_root_subtree_label != ROOT_SUBTREE_LABEL) {
            old_root_subtree_label = current_root_subtree_label;
            current_root_subtree_label = current_root_subtree_label >> subtree_depth;
            visited_subtree_labels.push_back(current_root_subtree_label);
            //cout << "Current remaining distance, root: " << current_remaining_distance << " " << current_root_subtree_label << "\n";
            for (unsigned long long partial_label = 0; partial_label < subtree_valence; ++partial_label) {
                neighboring_subtree_label = (current_root_subtree_label << subtree_depth) + partial_label;
                unsigned int recursive_remaining_distance = remaining_distance - label_distance(source_subtree_label, neighboring_subtree_label);
                if (neighboring_subtree_label != old_root_subtree_label)
                    add_descendant_subtree_labels_to_vector(
                            neighboring_subtree_label, tree_depth, subtree_depth,
                            max_depth, source_subtree_label,
                            recursive_remaining_distance, visited_subtree_labels);
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


vector< pair<unsigned long long, unsigned long long> > get_all_edges(
        unsigned long long *labels, unsigned int threshold, unsigned int tree_depth, unsigned int subtree_depth,
        vector<unsigned long long> subtree_labels, unordered_map<unsigned long long, unsigned int> map_subtree_label_to_index,
        unsigned long long *subtree_label_to_label_index_list, unsigned int *subtree_member_list_indices) {
    vector< pair<unsigned long long, unsigned long long> > edges;

    unsigned long long num_subtrees = subtree_labels.size();
    unsigned long long *first_subtree_member_labels, *second_subtree_member_labels;
    unsigned int num_first_subtree_member_labels, num_second_subtree_member_labels;
    vector<unsigned long long> visited_subtree_labels;
    vector<unsigned int> subtree_neighbors_accumulator;

    for (unsigned long long i = 0; i < num_subtrees; ++i) {
        //subtree_neighbors_accumulator.push_back(0);
        unsigned long long first_subtree_label = subtree_labels[i];
        first_subtree_member_labels = &subtree_label_to_label_index_list[subtree_member_list_indices[i]];
        num_first_subtree_member_labels = subtree_member_list_indices[i + 1] - subtree_member_list_indices[i];

        visited_subtree_labels.clear();
        int remaining_distance = threshold + 0 * subtree_depth;
        get_nearby_subtrees(first_subtree_label, tree_depth, subtree_depth, remaining_distance, visited_subtree_labels);
        //for (unsigned long long second_subtree_label : visited_subtree_labels) {
        for (unsigned long long j = 1; j < visited_subtree_labels.size(); ++j) {
            unsigned long long second_subtree_label = visited_subtree_labels[j];
            //++subtree_neighbors_accumulator[i];
            if (1 || second_subtree_label > first_subtree_label) {
                unsigned int second_subtree_index = map_subtree_label_to_index[second_subtree_label];
                second_subtree_member_labels = &subtree_label_to_label_index_list[subtree_member_list_indices[second_subtree_index]];
                num_second_subtree_member_labels = subtree_member_list_indices[second_subtree_index + 1];
                num_second_subtree_member_labels -= subtree_member_list_indices[second_subtree_index];

                for (unsigned int k = 0; k < num_first_subtree_member_labels; ++k) {
                    unsigned long long first_endpoint_label = first_subtree_member_labels[k];
                    for (unsigned int l = 0; l < num_second_subtree_member_labels; ++l) {
                        unsigned long long second_endpoint_label = second_subtree_member_labels[l];
                        if (label_distance(first_endpoint_label, second_endpoint_label) <= threshold)
                            edges.push_back(make_pair(first_endpoint_label, second_endpoint_label));
                    }
                }
            }
        }

        for (unsigned int k = 0; k < num_first_subtree_member_labels; ++k) {
            unsigned long long first_endpoint_label = first_subtree_member_labels[k];
            for (unsigned int l = k + 1; l < num_first_subtree_member_labels; ++l) {
                unsigned long long second_endpoint_label = first_subtree_member_labels[l];
                if (label_distance(first_endpoint_label, second_endpoint_label) <= threshold)
                    edges.push_back(pair<unsigned long long, unsigned long long>(first_endpoint_label, second_endpoint_label));
            }
        }
    }
    cout << "Average subtree neighborhood size: "
         << accumulate(subtree_neighbors_accumulator.begin(), subtree_neighbors_accumulator.end(), 0) * 1.0 / num_subtrees
         << endl << endl;

    return edges;
}
