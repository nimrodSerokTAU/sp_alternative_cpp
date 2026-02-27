#include "node.h"
#include <sstream>
#include <algorithm>

Node::Node(int node_id, const std::set<std::string>& keys, const std::vector<Node*>& children,
           double children_bl_sum, double branch_length)
    : id(node_id), keys(keys), father(nullptr), branch_length(branch_length),
      children(children), children_bl_sum(children_bl_sum),
      bl_sum_on_differentiator(0), bl_sum_on_non_differentiator(0),
      rank_from_root(-1), weight(0) {}

Node* Node::create_from_children(const std::vector<Node*>& children_list, int inx,
                                  std::vector<Node*>& storage) {
    std::set<std::string> combined_keys;
    double cbl_sum = 0;
    for (auto* child : children_list) {
        combined_keys.insert(child->keys.begin(), child->keys.end());
        cbl_sum += child->children_bl_sum + child->branch_length;
    }
    Node* node = new Node(inx, combined_keys, children_list, cbl_sum);
    storage.push_back(node);
    return node;
}

void Node::add_child_to_me(Node* child_node) {
    children.push_back(child_node);
    child_node->set_a_father(this);
}

void Node::set_a_father(Node* other_node) {
    father = other_node;
}

void Node::set_branch_length(double bl) {
    branch_length = bl;
}

std::string Node::get_keys_rooted_string() const {
    std::vector<std::string> sorted_keys(keys.begin(), keys.end());
    std::sort(sorted_keys.begin(), sorted_keys.end());
    std::string result;
    for (size_t i = 0; i < sorted_keys.size(); i++) {
        if (i > 0) result += ",";
        result += sorted_keys[i];
    }
    return result;
}

std::string Node::get_keys_unrooted_string(const std::set<std::string>& tree_keys,
                                            const std::string& differentiator_key) const {
    std::set<std::string> other_side;
    for (const auto& k : tree_keys) {
        if (keys.find(k) == keys.end()) {
            other_side.insert(k);
        }
    }

    const std::set<std::string>* use_set;
    if (keys.size() > other_side.size() ||
        (keys.size() == other_side.size() && other_side.count(differentiator_key))) {
        use_set = &other_side;
    } else {
        use_set = &keys;
    }

    if (use_set->size() < 2) {
        return "";
    }

    std::vector<std::string> sorted_keys(use_set->begin(), use_set->end());
    std::sort(sorted_keys.begin(), sorted_keys.end());
    std::string result;
    for (size_t i = 0; i < sorted_keys.size(); i++) {
        if (i > 0) result += ",";
        result += sorted_keys[i];
    }
    return result;
}

void Node::fill_newick() {
    if (children.empty()) {
        newick_part = *keys.begin() + ":" + std::to_string(branch_length);
    } else if (children.size() == 2) {
        newick_part = "(" + children[0]->newick_part + "," + children[1]->newick_part + "):" +
                      std::to_string(branch_length);
    }
}

void Node::set_parsimony_set(const std::set<std::string>& new_set) {
    parsimony_set = new_set;
}

std::vector<AdjEntry> Node::get_adj() const {
    std::vector<AdjEntry> res;
    for (auto* child : children) {
        res.push_back({child->id, child->branch_length});
    }
    if (father) {
        res.push_back({father->id, branch_length});
    }
    return res;
}

void Node::update_children_only(const std::vector<Node*>& children_list) {
    children = children_list;
}

void Node::set_rank_from_root(int rank) {
    rank_from_root = rank;
}

void Node::set_w_from_root(const std::vector<double>& w_list) {
    w_from_root = w_list;
}

void Node::update_data_from_children() {
    if (!children.empty()) {
        keys.clear();
        children_bl_sum = 0;
        for (auto* child : children) {
            keys.insert(child->keys.begin(), child->keys.end());
            children_bl_sum += child->children_bl_sum + child->branch_length;
        }
    }
}

void Node::set_weight_from_root() {
    double w = 0;
    std::vector<double> rev_w(w_from_root.rbegin(), w_from_root.rend());
    for (size_t index = 0; index < rev_w.size(); index++) {
        w += rev_w[index] / (index + 1);
    }
    weight = w;
}
