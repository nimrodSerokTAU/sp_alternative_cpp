#include "tree_stats.h"
#include "utils.h"
#include <algorithm>
#include <numeric>
#include <unordered_map>
#include <set>
#include <functional>

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

std::vector<int> calc_parsimony(const UnrootedTree& unrooted_tree,
    const std::vector<std::string>& aln,
    const std::vector<std::string>& names) {
    std::unordered_map<Node*, Node*> old_to_new;
    std::vector<std::unique_ptr<Node>> temp_nodes;
    temp_nodes.reserve(unrooted_tree.all_nodes.size());

    // Step 1: copy nodes
    for (auto& n : unrooted_tree.all_nodes) {
        auto copy = std::make_unique<Node>(n->id, n->keys, std::vector<Node*>{}, n->children_bl_sum, n->branch_length);
        old_to_new[n.get()] = copy.get();
        temp_nodes.push_back(std::move(copy));
    }

    // Step 2: fix children pointers
    for (auto& n : unrooted_tree.all_nodes) {
        Node* n_copy = old_to_new[n.get()];
        for (Node* c : n->children) {
            n_copy->children.push_back(old_to_new[c]);
        }
        if (n->father) n_copy->father = old_to_new[n->father];
    }

    // Step 3: get copied anchor
    Node* anchor_copy = old_to_new[unrooted_tree.anchor];

    // 3. If root has >2 children, create temporary binary root
    Node* temp_root = nullptr;
    if (anchor_copy->children.size() > 2) {
        auto new_node = std::make_unique<Node>(
            -1, std::set<std::string>{}, std::vector<Node*>{anchor_copy->children[0], anchor_copy->children[1]}, 0.0, 0.0);

        anchor_copy->children.erase(anchor_copy->children.begin(), anchor_copy->children.begin() + 2);
        anchor_copy->children.push_back(new_node.get());
        new_node->father = nullptr; // temporary root
        temp_nodes.push_back(std::move(new_node));
        temp_root = temp_nodes.back().get();
    }
    else {
        temp_root = anchor_copy;
    }

    // 4. Post-order traversal
    std::vector<Node*> postorder;
    std::function<void(Node*)> dfs = [&](Node* n) {
        for (Node* c : n->children) dfs(c);
        postorder.push_back(n);
        };
    dfs(temp_root);

    // 5. Map sequence names to indices
    std::unordered_map<std::string, int> seq_name_to_index;
    for (int i = 0; i < static_cast<int>(names.size()); i++)
        seq_name_to_index[names[i]] = i;

    // 6. Fitch algorithm
    std::vector<int> parsimony_per_col;
    int ncol = static_cast<int>(aln[0].size());
    for (int col = 0; col < ncol; col++) {
        std::unordered_map<Node*, std::set<char>> states;
        int score = 0;
        for (Node* node : postorder) {
            if (node->children.empty()) {
                char c = aln[seq_name_to_index[*node->keys.begin()]][col];
                states[node] = { c };
            }
            else {
                std::set<char> inter, uni;
                bool first = true;
                for (Node* c : node->children) {
                    if (first) { inter = states[c]; uni = states[c]; first = false; }
                    else {
                        std::set<char> temp;
                        std::set_intersection(inter.begin(), inter.end(),
                            states[c].begin(), states[c].end(),
                            std::inserter(temp, temp.begin()));
                        inter = temp;
                        uni.insert(states[c].begin(), states[c].end());
                    }
                }
                if (!inter.empty()) states[node] = inter;
                else { states[node] = uni; score++; }
            }
        }
        parsimony_per_col.push_back(score);
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
