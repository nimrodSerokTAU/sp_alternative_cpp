#include "rooted_tree.h"
#include <algorithm>
#include <cmath>
#include <deque>
#include <numeric>
#include <set>
#include <tuple>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include "utils.h"

RootedTree::RootedTree(Node* r,
    std::vector<std::unique_ptr<Node>> nodes,
    const std::set<std::string>& k)
{
    root = r;
    all_nodes = std::move(nodes);
    keys = k;
}


RootedTree::RootedTree(const UnrootedTree& unrooted, RootingMethods rooting_method) {
    RootingPoint rp;

    all_nodes.reserve(unrooted.all_nodes.size());
    keys = unrooted.keys;
    // copy nodes
    for (const auto& n : unrooted.all_nodes) {
        all_nodes.push_back(std::make_unique<Node>(*n));
    }

    // Find midpoint rooting location
    if (rooting_method == RootingMethods::LONGEST_PATH_MID) {
        rp = calc_mid_point(unrooted);
    }
    else {
        auto rooting_points = calc_min_differential_sum(unrooted);
        rp = find_shallowest_tree(unrooted, rooting_points);
    }

    Node* root_start = all_nodes[rp.start_id].get();
    Node* root_end = all_nodes[rp.end_id].get();



    // Create new root node
    std::set<std::string> root_keys = root_start->keys;
    root_keys.insert(root_end->keys.begin(), root_end->keys.end());

    auto new_root_ptr = std::make_unique<Node>(
        static_cast<int>(all_nodes.size()),
        root_keys,
        std::vector<int>{root_start->id, root_end->id},
        root_start->children_bl_sum + root_end->children_bl_sum,
        0.0
    );
    
    all_nodes.push_back(std::move(new_root_ptr));
    root = all_nodes.back().get();

    // Adjust branch lengths
    root_start->branch_length = rp.dist_from_start;
    root_end->branch_length = rp.dist_from_end;

    root_start->father_id = root->id;
    root_end->father_id = root->id;
    
	// Recalculate ranks and data for the affected subtree
    root->set_rank_from_root(0);
    std::vector<std::tuple<Node*, Node*, int>> nodes_to_recalc;
	std::vector<Node*> raw_nodes = get_raw_pointers_from_unique(all_nodes);


    nodes_to_recalc.push_back({ raw_nodes[rp.start_id], root, rp.end_id });
    nodes_to_recalc.push_back({ raw_nodes[rp.end_id], root, rp.start_id });

    while (!nodes_to_recalc.empty()) {
        auto [node, father, broke] = nodes_to_recalc.back();
        nodes_to_recalc.pop_back();
        recalc_tree_down(node, father, broke, nodes_to_recalc, raw_nodes);
    }

    // Sort by rank and update data from children (highest rank first)
    std::vector<Node*> sorted_by_rank = raw_nodes;
    std::sort(sorted_by_rank.begin(), sorted_by_rank.end(),
        [](const Node* a, const Node* b) { return a->rank_from_root > b->rank_from_root; });

    for (auto* node : sorted_by_rank) {
        node->update_data_from_children(raw_nodes);
    }
}

void RootedTree::calc_clustal_w(std::vector<Node*> raw_nodes) {
    std::deque<Node*> nodes_to_recalc;
    nodes_to_recalc.push_back(root);
    while (!nodes_to_recalc.empty()) {
        Node* node = nodes_to_recalc.front();
        nodes_to_recalc.pop_front();
        fill_nodes_w(node, nodes_to_recalc, raw_nodes);
    }
}

void RootedTree::calc_seq_w() {
	std::vector<Node*> raw_nodes = get_raw_pointers_from_unique(all_nodes);
    calc_clustal_w(raw_nodes);
    for (const auto& node : all_nodes) {
        if (node->children_ids.empty()) {
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
    for (const auto& adj : origin->get_adj(all_nodes)) {
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
            for (const auto& adj : next.node->get_adj(all_nodes)) {
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

std::vector<RootingPoint> calc_min_differential_sum(const UnrootedTree& unrooted) {
    
    std::vector<Node*> raw_nodes = get_raw_pointers_from_unique(unrooted.all_nodes);
    
    struct BranchInfo {
        int origin, dest;
        double bl, w_to_orig, w_to_dest, delta;
    };

    std::vector<double> all_bl = unrooted.get_branches_lengths_list();
    double total_branch_length = std::accumulate(all_bl.begin(), all_bl.end(), 0.0);
    double min_bl = *std::min_element(all_bl.begin(), all_bl.end());

    std::map<std::string, BranchInfo> all_branches;



    for (const auto& node : unrooted.all_nodes) {
        for (const auto& adj : node->get_adj(raw_nodes)) {
            int min_id = std::min(adj.node_id, node->id);
            int max_id = std::max(adj.node_id, node->id);
            std::string key = "o:" + std::to_string(min_id) + "-d:" + std::to_string(max_id);

            if (all_branches.find(key) == all_branches.end()) {
                double w_to_orig_val = sum_bl_up_to_node_id(raw_nodes[min_id], max_id, raw_nodes);
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
    for (int c_id : node->children_ids) {
        if (c_id == father->id) {
            father_is_child = true;
            break;
        }
    }
    if (father_is_child) {
        node->branch_length = father->branch_length;
    }

    auto adj = node->get_adj(all_nodes);
    std::vector<int> new_children;
    for (const auto& a : adj) {
        if (a.node_id != father->id && a.node_id != broke_id) {
            //for (auto* n : all_nodes) {
            //    if (n->id == a.node_id) {
                    //new_children.push_back(n->id);
            //        break;
            //    }
            //}
                    new_children.push_back(a.node_id);
        }
    }

    if (!new_children.empty()) {
        node->update_children_only(new_children);
        if (new_children.size() >= 1) {
            nodes_to_recalc.push_back({ all_nodes[new_children[0]], node, broke_id});
        }
        if (new_children.size() >= 2) {
            nodes_to_recalc.push_back({ all_nodes[new_children[1]], node, broke_id});
        }
    }
    node->father_id = father->id;
}

void fill_nodes_w(Node* node, std::deque<Node*>& nodes_to_recalc, std::vector<Node*>& all_nodes) {
    std::vector<double> w_list;
    if (node->father_id != -1) {
		Node* father = all_nodes[node->father_id];
        w_list = father->w_from_root;
    }
    w_list.push_back(node->branch_length);
    node->set_w_from_root(w_list);

    if (!node->children_ids.empty()) {
        nodes_to_recalc.push_back(all_nodes[node->children_ids[0]]);
        if (node->children_ids.size() > 1) {
            nodes_to_recalc.push_back(all_nodes[node->children_ids[1]]);
        }
    } else {
        node->set_weight_from_root();
    }
}


RootingPoint RootedTree::find_shallowest_tree(
    const UnrootedTree& ut,
    const std::vector<RootingPoint>& rooting_points)
{
    RootingPoint best = rooting_points[0];
    double best_depth = std::numeric_limits<double>::max();

    for (const auto& rp : rooting_points)
    {
        Node* start = ut.all_nodes[rp.start_id].get();
        Node* end = ut.all_nodes[rp.end_id].get();

        double depth = longest_dist_from_virtual_root(
            ut,
            start,
            end,
            rp.dist_from_start,
            rp.dist_from_end);

        if (depth < best_depth)
        {
            best_depth = depth;
            best = rp;
            best.longest_dist = depth;
        }
    }

    return best;
}

double longest_dist_from_virtual_root(
    const UnrootedTree& tree,
    Node* start,
    Node* end,
    double d_start,
    double d_end)
{
    std::deque<std::pair<Node*, double>> q;
    std::unordered_set<int> visited;

    double max_dist = 0.0;

    q.push_back({ start, d_start });
    q.push_back({ end,   d_end });

    visited.insert(start->id);
    visited.insert(end->id);

    std::vector<Node*> raw_nodes = get_raw_pointers_from_unique(tree.all_nodes);

    while (!q.empty())
    {
        auto [node, dist] = q.front();
        q.pop_front();

        max_dist = std::max(max_dist, dist);

        for (const auto& adj : node->get_adj(raw_nodes))
        {
            if (visited.count(adj.node_id))
                continue;

            Node* next = tree.all_nodes[adj.node_id].get();

            visited.insert(adj.node_id);
            q.push_back({ next, dist + adj.dist });
        }
    }

    return max_dist;
}
