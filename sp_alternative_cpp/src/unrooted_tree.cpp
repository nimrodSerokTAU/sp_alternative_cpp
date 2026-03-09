#include "unrooted_tree.h"
#include <fstream>
#include <algorithm>
#include <queue>
#include <stdexcept>

UnrootedTree::UnrootedTree(Node* anchor, const std::vector<Node*>& all_nodes)
    : anchor(anchor), all_nodes(all_nodes) {
    keys = anchor->keys;
    std::vector<std::string> sorted(keys.begin(), keys.end());
    std::sort(sorted.begin(), sorted.end());
    differentiator_key = sorted[0];
}

void UnrootedTree::take_ownership(Node* node) {
    owned_nodes.emplace_back(node);
}

void UnrootedTree::own_all(const std::vector<Node*>& nodes) {
    for (auto* n : nodes) {
        owned_nodes.emplace_back(n);
    }
}

UnrootedTree UnrootedTree::create_from_newick_file(const std::filesystem::path& path) {
    std::string newick_str = read_newick_from_file(path);
    return create_from_newick_str(newick_str);
}

UnrootedTree UnrootedTree::create_from_newick_str(const std::string& newick_str) {
    auto [root, all_nodes] = root_from_newick_str(newick_str);
    UnrootedTree tree(root, all_nodes);
    tree.own_all(all_nodes);
    return tree;
}

std::set<std::string> UnrootedTree::get_internal_edges_set() const {
    std::set<std::string> edges;
    for (auto* n : all_nodes) {
        if (!n->children.empty()) {
            std::string edge_str = n->get_keys_unrooted_string(keys, differentiator_key);
            if (!edge_str.empty()) {
                edges.insert(edge_str);
            }
        }
    }
    return edges;
}

int UnrootedTree::calc_rf(const UnrootedTree& other_tree) const {
    auto set_a = get_internal_edges_set();
    auto set_b = other_tree.get_internal_edges_set();
    // Symmetric difference
    int count = 0;
    for (const auto& e : set_a) {
        if (set_b.find(e) == set_b.end()) count++;
    }
    for (const auto& e : set_b) {
        if (set_a.find(e) == set_a.end()) count++;
    }
    return count;
}

std::vector<double> UnrootedTree::get_branches_lengths_list() const {
    std::vector<double> bl_list;
    for (auto* n : all_nodes) {
        if (n->father != nullptr) {
            bl_list.push_back(std::max(n->branch_length, TREE_EPSILON));
        }
    }
    return bl_list;
}

std::tuple<Node*, std::vector<PathEntry>, double> UnrootedTree::get_longest_path(Node* u) const {
    int nodes_count = static_cast<int>(all_nodes.size());
    std::vector<Node*> nodes_by_id(nodes_count + 1, nullptr);
    std::vector<double> distance(nodes_count + 1, -1.0);
    std::vector<std::vector<PathEntry>> path(nodes_count + 1);

    distance[u->id] = 0;
    nodes_by_id[u->id] = u;

    std::deque<Node*> queue_nodes;
    queue_nodes.push_back(u);

    while (!queue_nodes.empty()) {
        Node* front = queue_nodes.front();
        queue_nodes.pop_front();

        for (const auto& adj : front->get_adj()) {
            if (adj.node_id >= 0 && adj.node_id <= nodes_count && nodes_by_id[adj.node_id] == nullptr) {
                // Find the actual node
                Node* adj_node = nullptr;
                for (auto* n : all_nodes) {
                    if (n->id == adj.node_id) {
                        adj_node = n;
                        break;
                    }
                }
                if (!adj_node) continue;

                nodes_by_id[adj.node_id] = adj_node;
                distance[adj.node_id] = distance[front->id] + adj.dist;
                path[adj.node_id] = path[front->id];
                path[adj.node_id].push_back({front->id, adj.node_id, adj.dist});
                queue_nodes.push_back(adj_node);
            }
        }
    }

    double max_dist = 0;
    int node_index = -1;
    for (int i = 0; i < nodes_count; i++) {
        if (distance[i] > max_dist) {
            max_dist = distance[i];
            node_index = i;
        }
    }

    return {nodes_by_id[node_index], path[node_index], max_dist};
}

std::pair<std::vector<PathEntry>, double> UnrootedTree::longest_path() const {
    auto [node1, path1, dist1] = get_longest_path(all_nodes[0]);
    auto [node2, path2, dist2] = get_longest_path(node1);
    return {path2, dist2};
}

double UnrootedTree::get_longest_dist_to(Node* dest) const {
    int nodes_count = static_cast<int>(all_nodes.size());
    std::vector<Node*> nodes_by_id(nodes_count + 1, nullptr);
    std::vector<double> distance(nodes_count + 1, -1.0);

    distance[dest->id] = 0;
    nodes_by_id[dest->id] = dest;

    std::deque<Node*> queue_nodes;
    queue_nodes.push_back(dest);

    while (!queue_nodes.empty()) {
        Node* front = queue_nodes.front();
        queue_nodes.pop_front();

        for (const auto& adj : front->get_adj()) {
            if (adj.node_id >= 0 && adj.node_id <= nodes_count && nodes_by_id[adj.node_id] == nullptr) {
                Node* adj_node = nullptr;
                for (auto* n : all_nodes) {
                    if (n->id == adj.node_id) {
                        adj_node = n;
                        break;
                    }
                }
                if (!adj_node) continue;

                nodes_by_id[adj.node_id] = adj_node;
                distance[adj.node_id] = distance[front->id] + adj.dist;
                queue_nodes.push_back(adj_node);
            }
        }
    }

    return *std::max_element(distance.begin(), distance.end());
}

// Free functions

std::string read_newick_from_file(const std::filesystem::path& input_file_path) {
    std::ifstream in_file(input_file_path);
    if (!in_file.is_open()) {
        throw std::runtime_error("Cannot open newick file: " + input_file_path.string());
    }
    std::string line;
    if (std::getline(in_file, line)) {
        while (!line.empty() && (line.back() == '\r' || line.back() == '\n' || line.back() == ' ')) {
            line.pop_back();
        }
        return line;
    }
    return "";
}

Node* create_node_from_children(std::vector<std::vector<Node*>>& open_nodes_per_level,
                                 int level, double branch_length, int node_inx,
                                 std::vector<Node*>& storage) {
    int child_level = level + 1;
    if (child_level >= static_cast<int>(open_nodes_per_level.size())) {
        return nullptr;
    }

    std::set<std::string> node_keys;
    for (auto* child : open_nodes_per_level[child_level]) {
        node_keys.insert(child->keys.begin(), child->keys.end());
    }

    std::vector<Node*> child_list = open_nodes_per_level[child_level];
    Node* current_node = new Node(node_inx, node_keys, child_list, 0, branch_length);
    storage.push_back(current_node);

    for (auto* child : open_nodes_per_level[child_level]) {
        child->set_a_father(current_node);
    }

    return current_node;
}

std::pair<Node*, std::vector<Node*>> create_a_tree_from_newick(const std::string& newick) {
    std::vector<Node*> all_nodes;
    // Use a large enough vector for levels
    std::vector<std::vector<Node*>> open_nodes_per_level(100);
    int level = 0;
    std::string current_key;
    std::string branch_length_str;
    size_t i = 0;

    while (i < newick.size()) {
        if (newick[i] == ' ') {
            i++;
        } else if (newick[i] == '(') {
            level++;
            if (level >= static_cast<int>(open_nodes_per_level.size())) {
                open_nodes_per_level.resize(level + 10);
            }
            open_nodes_per_level[level].clear();
            i++;
        } else if (newick[i] == ':') {
            branch_length_str.clear();
            i++;
            while (i < newick.size() && newick[i] != ')' && newick[i] != ',') {
                branch_length_str += newick[i];
                i++;
            }
        } else if (newick[i] == ',' || newick[i] == ')') {
            double bl = branch_length_str.empty() ? 0.0 : std::stod(branch_length_str);
            if (!current_key.empty()) {
                Node* current_node = new Node(static_cast<int>(all_nodes.size()),
                                               {current_key}, {}, 0, bl);
                all_nodes.push_back(current_node);
                open_nodes_per_level[level].push_back(current_node);
                current_key.clear();
            } else {
                Node* current_node = create_node_from_children(open_nodes_per_level, level,
                                                                bl, static_cast<int>(all_nodes.size()),
                                                                all_nodes);
                open_nodes_per_level[level + 1].clear();
                open_nodes_per_level[level].push_back(current_node);
            }
            if (newick[i] == ')') {
                level--;
            }
            i++;
        } else if (newick[i] == ';') {
            double bl = branch_length_str.empty() ? 0.0 : std::stod(branch_length_str);
            Node* current_node = create_node_from_children(open_nodes_per_level, 0,
                                                            bl, static_cast<int>(all_nodes.size()),
                                                            all_nodes);
            return {current_node, all_nodes};
        } else {
            current_key.clear();
            while (i < newick.size() && newick[i] != ')' && newick[i] != ',' && newick[i] != ':') {
                current_key += newick[i];
                i++;
            }
        }
    }

    return {nullptr, all_nodes};
}

std::pair<Node*, std::vector<Node*>> root_from_newick_str(const std::string& newick_str) {
    auto [root, all_nodes] = create_a_tree_from_newick(newick_str);

    if (root->children.size() == 3) {
        return {root, all_nodes};
    }

    if (root->children.size() == 2) {
        // Sort children: by branch_length ascending, then by number of children descending
        std::sort(root->children.begin(), root->children.end(),
                  [](const Node* a, const Node* b) {
                      if (a->children.size() != b->children.size())
                          return a->children.size() > b->children.size();
                      return a->branch_length < b->branch_length;
                  });

        if (root->children[0]->children.size() == 2) {
            std::vector<Node*> res;
            res.push_back(root->children[0]->children[0]);
            res.push_back(root->children[0]->children[1]);
            res.push_back(root->children[1]);

            std::set<std::string> combined_keys;
            double cbl_sum = 0;
            for (auto* child : res) {
                combined_keys.insert(child->keys.begin(), child->keys.end());
                cbl_sum += child->children_bl_sum + child->branch_length;
            }
            Node* anchor = new Node(static_cast<int>(all_nodes.size()), combined_keys, res, cbl_sum);
            all_nodes.push_back(anchor);
            return {anchor, all_nodes};
        }
    }

    return {root, all_nodes};
}
