#include "neighbor_joining.h"
#include <limits>
#include <numeric>
#include <algorithm>

using namespace std;

NeighborJoining::NeighborJoining(
    const std::vector<std::vector<double>>& distanceMatrix,
    std::vector<std::unique_ptr<Node>> initial_nodes)
{
    distance_matrix = distanceMatrix;

    // move ownership of nodes
    all_nodes = std::move(initial_nodes);

    // build working_nodes from the owned nodes
    working_nodes.reserve(all_nodes.size());
    for (auto& node : all_nodes) {
        working_nodes.push_back(node.get());
    }

    // initialize q_matrix size
    q_matrix.resize(distance_matrix.size(),
        std::vector<double>(distance_matrix.size(), 0.0));

    tree_res = build_tree();
}

std::vector<std::vector<double>> NeighborJoining::calc_q_matrix() const {
    int n = static_cast<int>(working_nodes.size());
    std::vector<std::vector<double>> qm(distance_matrix.size() - 1, std::vector<double>(n, 0.0));

    for (int i = 0; i < static_cast<int>(distance_matrix.size()) - 1; i++) {
        double sum_i = std::accumulate(distance_matrix[i].begin(), distance_matrix[i].end(), 0.0);
        for (int j = i + 1; j < static_cast<int>(distance_matrix.size()); j++) {
            double sum_j = std::accumulate(distance_matrix[j].begin(), distance_matrix[j].end(), 0.0);
            qm[i][j] = (n - 2) * distance_matrix[i][j] - sum_i - sum_j;
        }
    }
    return qm;
}

NeighborJoining::ClosestPair NeighborJoining::find_closest_pair() const {
    int f_inx = -1, s_inx = -1;
    double min_dist = std::numeric_limits<double>::max();

    for (int i = 0; i < static_cast<int>(q_matrix.size()); i++) {
        for (int j = i + 1; j < static_cast<int>(q_matrix[i].size()); j++) {
            if (q_matrix[i][j] < min_dist) {
                min_dist = q_matrix[i][j];
                f_inx = i;
                s_inx = j;
            }
        }
    }
    return {f_inx, s_inx, min_dist};
}

std::pair<double, double> NeighborJoining::find_delta(int f_inx, int s_inx) const {
    int n = static_cast<int>(working_nodes.size());
    double sum_f = std::accumulate(distance_matrix[f_inx].begin(), distance_matrix[f_inx].end(), 0.0);
    double sum_s = std::accumulate(distance_matrix[s_inx].begin(), distance_matrix[s_inx].end(), 0.0);

    double delta_f = 0.5 * distance_matrix[f_inx][s_inx] +
                     (sum_f - sum_s) / (2.0 * (n - 2));
    double delta_s = distance_matrix[f_inx][s_inx] - delta_f;
    return {delta_f, delta_s};
}

void NeighborJoining::merge_two_clusters() {
    auto [f_inx, s_inx, min_dist] = find_closest_pair();
    auto [delta_f, delta_s] = find_delta(f_inx, s_inx);

    std::vector<std::vector<double>> matrix;
    for (int r = 0; r < static_cast<int>(distance_matrix.size()); r++) {
        std::vector<double> row = distance_matrix[r];
        if (r == f_inx) {
            row[f_inx] = 0;
        } else {
            row[f_inx] = (row[f_inx] + row[s_inx] - distance_matrix[f_inx][s_inx]) / 2.0;
        }
        row.erase(row.begin() + s_inx);

        if (r != s_inx) {
            matrix.push_back(row);
        } else {
            // Merge into f_inx row
            for (int n_idx = 0; n_idx < static_cast<int>(row.size()); n_idx++) {
                if (n_idx != f_inx) {
                    matrix[f_inx][n_idx] = (matrix[f_inx][n_idx] + row[n_idx] -
                                             distance_matrix[f_inx][s_inx]) / 2.0;
                }
            }
        }
    }

    distance_matrix = matrix;

    working_nodes[f_inx]->set_branch_length(delta_f);
    working_nodes[s_inx]->set_branch_length(delta_s);

    // Create new node
    set<string> keysOfNewNode = working_nodes[f_inx]->keys;
    keysOfNewNode.insert(working_nodes[s_inx]->keys.begin(), working_nodes[s_inx]->keys.end());
	double children_bl_sum = delta_f + delta_s;

    auto new_node = std::make_unique<Node>(static_cast<int>(all_nodes.size()),
                                            keysOfNewNode,
                                            std::vector<Node*>{working_nodes[f_inx], working_nodes[s_inx]},
                                            children_bl_sum,
                                            0);

    Node* new_node_ptr = new_node.get();
    all_nodes.push_back(std::move(new_node));
    working_nodes[f_inx]->set_a_father(new_node_ptr);
    working_nodes[s_inx]->set_a_father(new_node_ptr);

    
    working_nodes[f_inx] = new_node_ptr;
    working_nodes.erase(working_nodes.begin() + s_inx);
}

Node* NeighborJoining::merge_last_three() {
    auto [delta_f, delta_s] = find_delta(0, 1);
    double delta_t = distance_matrix[0][2] - delta_f;
    double deltas[3] = {delta_f, delta_s, delta_t};

    for (int i = 0; i < 3; i++) {
        working_nodes[i]->set_branch_length(deltas[i]);
    }

    auto new_node = std::make_unique<Node>(static_cast<int>(all_nodes.size()),
                                            std::set<std::string>{},
                                            std::vector<Node*>{working_nodes[0], working_nodes[1], working_nodes[2]}, 0);
    for (auto* child : new_node->children) {
        new_node->keys.insert(child->keys.begin(), child->keys.end());
        new_node->children_bl_sum += child->children_bl_sum + child->branch_length;
    }

    for (int i = 0; i < 3; i++) {
        working_nodes[i]->set_a_father(new_node.get());
    }
    all_nodes.push_back(std::move(new_node));
    return new_node.get();
}

UnrootedTree NeighborJoining::build_tree() {
    while (working_nodes.size() > 3) {
        q_matrix = calc_q_matrix();
        merge_two_clusters();
    }
    q_matrix = calc_q_matrix();
    Node* root = merge_last_three();
    return UnrootedTree(root, std::move(all_nodes));
}
