#include "rooted_tree.h"
#include <algorithm>
#include <cmath>
#include <deque>
#include <numeric>
#include <set>
#include <tuple>
#include <map>
#include <unordered_map>

// Deep copy a vector of nodes preserving parent/child relationships by id
std::vector<Node*> deep_copy_nodes(const std::vector<Node*>& src,
                                    std::vector<std::unique_ptr<Node>>& storage) {
    std::unordered_map<int, Node*> id_map;
    std::vector<Node*> result;

    // First pass: create copies
    for (auto* n : src) {
        auto copy = std::make_unique<Node>(n->id, n->keys, std::vector<Node*>{},
                                            n->children_bl_sum, n->branch_length);
        copy->rank_from_root = n->rank_from_root;
        copy->w_from_root = n->w_from_root;
        copy->weight = n->weight;
        copy->parsimony_set = n->parsimony_set;
        id_map[n->id] = copy.get();
        result.push_back(copy.get());
        storage.push_back(std::move(copy));
    }

    // Second pass: fix relationships
    for (size_t i = 0; i < src.size(); i++) {
        for (auto* child : src[i]->children) {
            auto it = id_map.find(child->id);
            if (it != id_map.end()) {
                result[i]->children.push_back(it->second);
            }
        }
        if (src[i]->father) {
            auto it = id_map.find(src[i]->father->id);
            if (it != id_map.end()) {
                result[i]->father = it->second;
            }
        }
    }

    return result;
}

RootedTree::RootedTree(Node* root, const std::vector<Node*>& all_nodes,
                        const std::set<std::string>& keys)
    : root(root), all_nodes(all_nodes), keys(keys) {}

RootedTree RootedTree::root_tree(const UnrootedTree& unrooted, RootingMethods rooting_method) {
    std::vector<std::unique_ptr<Node>> storage;
    std::vector<Node*> all_nodes = deep_copy_nodes(unrooted.all_nodes, storage);
    std::set<std::string> keys = unrooted.keys;

    RootingPoint rooting_point;
    if (rooting_method == RootingMethods::LONGEST_PATH_MID) {
        // Need to use the copy for midpoint calculation
        // Create a temporary unrooted tree from copied nodes
        // Find the anchor equivalent in copied nodes
        Node* copied_anchor = nullptr;
        for (auto* n : all_nodes) {
            if (n->id == unrooted.anchor->id) {
                copied_anchor = n;
                break;
            }
        }
        UnrootedTree temp_tree(copied_anchor, all_nodes);
        rooting_point = calc_mid_point(temp_tree);
    } else {
        Node* copied_anchor = nullptr;
        for (auto* n : all_nodes) {
            if (n->id == unrooted.anchor->id) {
                copied_anchor = n;
                break;
            }
        }
        UnrootedTree temp_tree(copied_anchor, all_nodes);
        auto rooting_points = calc_min_differential_sum(temp_tree, all_nodes);
        rooting_point = find_shallowest_tree(all_nodes, rooting_points);
    }

    auto [new_root, updated_nodes] = create_root(all_nodes, rooting_point, storage);

    RootedTree tree(new_root, updated_nodes, keys);
    tree.owned_nodes = std::move(storage);
    return tree;
}

void RootedTree::calc_clustal_w() {
    std::deque<Node*> nodes_to_recalc;
    nodes_to_recalc.push_back(root);
    while (!nodes_to_recalc.empty()) {
        Node* node = nodes_to_recalc.front();
        nodes_to_recalc.pop_front();
        fill_nodes_w(node, nodes_to_recalc);
    }
}

void RootedTree::calc_seq_w() {
    calc_clustal_w();
    for (auto* node : all_nodes) {
        if (node->children.empty()) {
            seq_weight_dict[*node->keys.begin()] = node->weight;
        }
    }
}

RootingPoint calc_mid_point(const UnrootedTree& unrooted) {
    auto [path, max_dist] = unrooted.longest_path();
    double half_length = 0;
    for (const auto& b : path) {
        if (half_length + b.dist > max_dist / 2.0) {
            double dist_from_start = max_dist / 2.0 - half_length;
            return {b.start_id, b.end_id, dist_from_start, b.dist - dist_from_start, 0};
        }
        half_length += b.dist;
    }
    // Fallback: use last edge
    if (!path.empty()) {
        return {path.back().start_id, path.back().end_id, path.back().dist / 2.0, path.back().dist / 2.0, 0};
    }
    return {0, 0, 0, 0, 0};
}

double sum_bl_up_to_node_id(Node* origin, int dest_id, const std::vector<Node*>& all_nodes) {
    double w_to_orig = 0;
    std::set<int> visited_node_ids = {dest_id, origin->id};

    struct QueueEntry { Node* node; double dist; };
    std::vector<QueueEntry> queue;
    for (const auto& adj : origin->get_adj()) {
        if (visited_node_ids.find(adj.node_id) == visited_node_ids.end()) {
            // Find actual node
            Node* adj_node = nullptr;
            for (auto* n : all_nodes) {
                if (n->id == adj.node_id) { adj_node = n; break; }
            }
            if (adj_node) queue.push_back({adj_node, adj.dist});
        }
    }

    while (!queue.empty()) {
        auto next = queue.back();
        queue.pop_back();
        if (visited_node_ids.find(next.node->id) == visited_node_ids.end()) {
            w_to_orig += next.dist;
            visited_node_ids.insert(next.node->id);
            for (const auto& adj : next.node->get_adj()) {
                if (visited_node_ids.find(adj.node_id) == visited_node_ids.end()) {
                    Node* adj_node = nullptr;
                    for (auto* n : all_nodes) {
                        if (n->id == adj.node_id) { adj_node = n; break; }
                    }
                    if (adj_node) queue.push_back({adj_node, adj.dist});
                }
            }
        }
    }
    return w_to_orig;
}

RootingPoint calc_potential_root_on_branch(double bl, double w_to_orig, double w_to_dest,
                                            double delta, int origin, int dest, double min_bl) {
    min_bl /= 10.0;
    double dist_from_start;
    if (delta > 0) {
        dist_from_start = bl / 2.0 - (w_to_orig - w_to_dest) / 2.0;
    } else if (w_to_orig > w_to_dest) {
        dist_from_start = min_bl;
    } else {
        dist_from_start = bl - min_bl;
    }
    return {origin, dest, dist_from_start, bl - dist_from_start, 0};
}

std::vector<RootingPoint> calc_min_differential_sum(const UnrootedTree& unrooted,
                                                      std::vector<Node*>& all_nodes) {
    struct BranchInfo {
        int origin, dest;
        double bl, w_to_orig, w_to_dest, delta;
    };

    std::vector<double> all_bl = unrooted.get_branches_lengths_list();
    double total_branch_length = std::accumulate(all_bl.begin(), all_bl.end(), 0.0);
    double min_bl = *std::min_element(all_bl.begin(), all_bl.end());

    std::map<std::string, BranchInfo> all_branches;

    for (auto* node : all_nodes) {
        for (const auto& adj : node->get_adj()) {
            int min_id = std::min(adj.node_id, node->id);
            int max_id = std::max(adj.node_id, node->id);
            std::string key = "o:" + std::to_string(min_id) + "-d:" + std::to_string(max_id);

            if (all_branches.find(key) == all_branches.end()) {
                double w_to_orig_val = sum_bl_up_to_node_id(all_nodes[min_id], max_id, all_nodes);
                double w_to_dest_val = total_branch_length - w_to_orig_val - adj.dist;
                double delta_val = adj.dist - std::abs(w_to_dest_val - w_to_orig_val);

                all_branches[key] = {min_id, max_id, adj.dist, w_to_orig_val, w_to_dest_val, delta_val};
            }
        }
    }

    std::vector<RootingPoint> res;
    for (const auto& [key, b] : all_branches) {
        if (b.delta > 0) {
            res.push_back(calc_potential_root_on_branch(b.bl, b.w_to_orig, b.w_to_dest,
                                                         b.delta, b.origin, b.dest, min_bl));
        }
    }

    if (res.empty()) {
        // Find branch with max delta
        auto best = std::max_element(all_branches.begin(), all_branches.end(),
                                      [](const auto& a, const auto& b) { return a.second.delta < b.second.delta; });
        res.push_back(calc_potential_root_on_branch(best->second.bl, best->second.w_to_orig,
                                                     best->second.w_to_dest, best->second.delta,
                                                     best->second.origin, best->second.dest, min_bl));
    }

    return res;
}

void recalc_tree_down(Node* node, Node* father, int broke_id,
                       std::vector<std::tuple<Node*, Node*, int>>& nodes_to_recalc,
                       std::vector<Node*>& all_nodes) {
    node->set_rank_from_root(father->rank_from_root + 1);

    // Check if father is in node's children
    bool father_is_child = false;
    for (auto* c : node->children) {
        if (c->id == father->id) {
            father_is_child = true;
            break;
        }
    }
    if (father_is_child) {
        node->branch_length = father->branch_length;
    }

    auto adj = node->get_adj();
    std::vector<Node*> new_children;
    for (const auto& a : adj) {
        if (a.node_id != father->id && a.node_id != broke_id) {
            for (auto* n : all_nodes) {
                if (n->id == a.node_id) {
                    new_children.push_back(n);
                    break;
                }
            }
        }
    }

    if (!new_children.empty()) {
        node->update_children_only(new_children);
        if (new_children.size() >= 1) {
            nodes_to_recalc.push_back({new_children[0], node, broke_id});
        }
        if (new_children.size() >= 2) {
            nodes_to_recalc.push_back({new_children[1], node, broke_id});
        }
    }
    node->father = father;
}

void fill_nodes_w(Node* node, std::deque<Node*>& nodes_to_recalc) {
    std::vector<double> w_list;
    if (node->father != nullptr) {
        w_list = node->father->w_from_root;
    }
    w_list.push_back(node->branch_length);
    node->set_w_from_root(w_list);

    if (!node->children.empty()) {
        nodes_to_recalc.push_back(node->children[0]);
        if (node->children.size() > 1) {
            nodes_to_recalc.push_back(node->children[1]);
        }
    } else {
        node->set_weight_from_root();
    }
}

std::pair<Node*, std::vector<Node*>> create_root(std::vector<Node*>& all_nodes,
                                                   const RootingPoint& rooting_point,
                                                   std::vector<std::unique_ptr<Node>>& storage) {
    int new_root_id = static_cast<int>(all_nodes.size());
    int start_id = rooting_point.start_id;
    int end_id = rooting_point.end_id;

    auto new_root_ptr = std::make_unique<Node>(new_root_id, std::set<std::string>{},
                                                std::vector<Node*>{all_nodes[start_id], all_nodes[end_id]},
                                                0);
    Node* new_root = new_root_ptr.get();
    new_root->set_rank_from_root(0);
    storage.push_back(std::move(new_root_ptr));
    all_nodes.push_back(new_root);

    std::vector<std::tuple<Node*, Node*, int>> nodes_to_recalc;
    nodes_to_recalc.push_back({all_nodes[start_id], new_root, end_id});
    nodes_to_recalc.push_back({all_nodes[end_id], new_root, start_id});

    while (!nodes_to_recalc.empty()) {
        auto [node, father, broke] = nodes_to_recalc.back();
        nodes_to_recalc.pop_back();
        recalc_tree_down(node, father, broke, nodes_to_recalc, all_nodes);
    }

    all_nodes[start_id]->set_branch_length(rooting_point.dist_from_start);
    all_nodes[end_id]->set_branch_length(rooting_point.dist_from_end);

    // Sort by rank and update data from children (highest rank first)
    std::vector<Node*> sorted_by_rank = all_nodes;
    std::sort(sorted_by_rank.begin(), sorted_by_rank.end(),
              [](const Node* a, const Node* b) { return a->rank_from_root > b->rank_from_root; });

    for (auto* node : sorted_by_rank) {
        node->update_data_from_children();
    }

    return {new_root, all_nodes};
}

RootingPoint find_shallowest_tree(std::vector<Node*>& all_nodes,
                                   std::vector<RootingPoint>& rooting_points) {
    if (rooting_points.size() == 1) {
        return rooting_points[0];
    }

    for (auto& rp : rooting_points) {
        std::vector<std::unique_ptr<Node>> temp_storage;
        auto my_nodes = deep_copy_nodes(all_nodes, temp_storage);

        // Add a node between
        int start_id = rp.start_id;
        int end_id = rp.end_id;

        bool start_father_is_end = (my_nodes[start_id]->father &&
                                     my_nodes[start_id]->father->id == end_id);

        Node *child, *parent;
        double dist_child, dist_parent;
        if (start_father_is_end) {
            child = my_nodes[start_id];
            parent = my_nodes[end_id];
            dist_child = rp.dist_from_start;
            dist_parent = rp.dist_from_end;
        } else {
            child = my_nodes[end_id];
            parent = my_nodes[start_id];
            dist_child = rp.dist_from_end;
            dist_parent = rp.dist_from_start;
        }

        // Create new node between child and parent
        std::set<std::string> new_keys = child->keys;
        std::vector<Node*> child_vec = {child};
        auto* new_node = new Node(static_cast<int>(my_nodes.size()),
                                   new_keys, child_vec, 0);
        temp_storage.emplace_back(new_node);
        new_node->set_a_father(parent);
        new_node->set_branch_length(dist_parent);

        // Update parent's children
        std::vector<Node*> new_parent_children;
        for (auto* c : parent->children) {
            if (c->id != child->id) new_parent_children.push_back(c);
        }
        new_parent_children.push_back(new_node);
        parent->children = new_parent_children;

        child->set_a_father(new_node);
        child->set_branch_length(dist_child);
        my_nodes.push_back(new_node);

        // Create a sentinel node for root measurement
        std::set<std::string> empty_keys;
        std::vector<Node*> empty_children;
        auto* sentinel = new Node(static_cast<int>(my_nodes.size()),
                                   empty_keys, empty_children, 0);
        temp_storage.emplace_back(sentinel);
        my_nodes.push_back(sentinel);

        // Find anchor equivalent
        Node* anchor = my_nodes.size() >= 3 ? my_nodes[my_nodes.size() - 3] : my_nodes[0];
        UnrootedTree temp_tree(anchor, my_nodes);
        rp.longest_dist = temp_tree.get_longest_dist_to(my_nodes[my_nodes.size() - 2]);
    }

    std::sort(rooting_points.begin(), rooting_points.end(),
              [](const RootingPoint& a, const RootingPoint& b) {
                  return a.longest_dist < b.longest_dist;
              });

    return rooting_points[0];
}
