#include "tree_stats.h"
#include "utils.h"
#include <algorithm>
#include <numeric>
#include <unordered_map>
#include <set>

TreeStats::TreeStats(const std::string& code, int taxa_num, int msa_length)
    : BasicStats(code, taxa_num, msa_length,
                  {"code", "bl_sum", "bl_mean", "bl_25_pct", "bl_75_pct", "bl_max", "bl_min",
                   "parsimony_mean", "parsimony_sum", "parsimony_max", "parsimony_min",
                   "parsimony_25_pct", "parsimony_75_pct"}),
      bl_mean(-1), bl_sum(-1), bl_max(-1), bl_min(-1), bl_25_pct(-1), bl_75_pct(-1),
      parsimony_mean(-1), parsimony_sum(-1), parsimony_max(-1), parsimony_min(-1),
      parsimony_25_pct(-1), parsimony_75_pct(-1) {}

void TreeStats::set_tree_stats(const std::vector<double>& bl_list, const UnrootedTree& tree,
                                 const std::vector<std::string>& aln, const std::vector<std::string>& names) {
    std::vector<double> bl = bl_list;
    bl_sum = std::accumulate(bl.begin(), bl.end(), 0.0);
    bl_mean = bl_sum / bl.size();
    bl_max = *std::max_element(bl.begin(), bl.end());
    bl_min = *std::min_element(bl.begin(), bl.end());
    bl_25_pct = calc_percentile(bl, 25);
    bl_75_pct = calc_percentile(bl, 75);

    auto parsimony_list = calc_parsimony(tree, aln, names);
    std::vector<double> pars_double(parsimony_list.begin(), parsimony_list.end());
    parsimony_sum = std::accumulate(pars_double.begin(), pars_double.end(), 0.0);
    parsimony_mean = parsimony_sum / pars_double.size();
    parsimony_max = *std::max_element(pars_double.begin(), pars_double.end());
    parsimony_min = *std::min_element(pars_double.begin(), pars_double.end());
    parsimony_25_pct = calc_percentile(pars_double, 25);
    parsimony_75_pct = calc_percentile(pars_double, 75);
}

std::vector<int> calc_parsimony(const UnrootedTree& unrooted_tree, const std::vector<std::string>& aln,
                                 const std::vector<std::string>& names) {
    std::vector<int> parsimony_per_col;

    // Create temporary nodes for rooting
    // new_node = children[0], children[1]
    // new_root = new_node, children[2]
    Node new_node(-1, {}, {unrooted_tree.anchor->children[0], unrooted_tree.anchor->children[1]}, 0);
    // Fill keys
    for (auto* c : new_node.children) {
        new_node.keys.insert(c->keys.begin(), c->keys.end());
    }
    Node new_root(-2, {}, {&new_node, unrooted_tree.anchor->children[2]}, 0);
    for (auto* c : new_root.children) {
        new_root.keys.insert(c->keys.begin(), c->keys.end());
    }

    // Build nodes_order: all_nodes except last (anchor), sorted by id, plus new_node, new_root
    std::vector<Node*> nodes_order;
    for (int i = 0; i < static_cast<int>(unrooted_tree.all_nodes.size()) - 1; i++) {
        nodes_order.push_back(unrooted_tree.all_nodes[i].get());
    }
    std::sort(nodes_order.begin(), nodes_order.end(), [](Node* a, Node* b) { return a->id < b->id; });
    nodes_order.push_back(&new_node);
    nodes_order.push_back(&new_root);

    std::unordered_map<std::string, int> seq_name_to_index;
    for (int i = 0; i < static_cast<int>(names.size()); i++) {
        seq_name_to_index[names[i]] = i;
    }

    int seq_len = static_cast<int>(aln[0].size());
    for (int col_index = 0; col_index < seq_len; col_index++) {
        int col_counter = 0;
        for (auto* n : nodes_order) {
            if (n->children.empty()) {
                std::string key = *n->keys.begin();
                int seq_index = seq_name_to_index[key];
                std::string ch(1, aln[seq_index][col_index]);
                n->set_parsimony_set({ch});
            } else {
                const auto& set_a = n->children[0]->parsimony_set;
                const auto& set_b = n->children[1]->parsimony_set;
                std::set<std::string> intersection;
                std::set_intersection(set_a.begin(), set_a.end(), set_b.begin(), set_b.end(),
                                      std::inserter(intersection, intersection.begin()));
                if (!intersection.empty()) {
                    n->set_parsimony_set(intersection);
                } else {
                    std::set<std::string> union_set;
                    std::set_union(set_a.begin(), set_a.end(), set_b.begin(), set_b.end(),
                                   std::inserter(union_set, union_set.begin()));
                    n->set_parsimony_set(union_set);
                    col_counter++;
                }
            }
        }
        parsimony_per_col.push_back(col_counter);
    }
    return parsimony_per_col;
}

std::vector<StatValue> TreeStats::get_my_features_as_list() const {
    return {
        StatValue(code),
        StatValue(bl_sum), StatValue(bl_mean), StatValue(bl_25_pct), StatValue(bl_75_pct),
        StatValue(bl_max), StatValue(bl_min),
        StatValue(parsimony_mean), StatValue(parsimony_sum), StatValue(parsimony_max), StatValue(parsimony_min),
        StatValue(parsimony_25_pct), StatValue(parsimony_75_pct)
    };
}
