#pragma once
#include "node.h"
#include <string>
#include <vector>
#include <set>
#include <memory>
#include <tuple>
#include <filesystem>

struct PathEntry {
    int start_id;
    int end_id;
    double dist;
};

class UnrootedTree {
public:
    Node* anchor;
    std::vector<Node*> all_nodes;
    std::set<std::string> keys;
    std::string differentiator_key;

    // Memory management: owns all dynamically allocated nodes
    std::vector<std::unique_ptr<Node>> owned_nodes;

    UnrootedTree(Node* anchor, const std::vector<Node*>& all_nodes);

    static UnrootedTree create_from_newick_file(const std::filesystem::path& path);
    static UnrootedTree create_from_newick_str(const std::string& newick_str);

    std::set<std::string> get_internal_edges_set() const;
    int calc_rf(const UnrootedTree& other_tree) const;
    std::vector<double> get_branches_lengths_list() const;

    std::tuple<Node*, std::vector<PathEntry>, double> get_longest_path(Node* u) const;
    std::pair<std::vector<PathEntry>, double> longest_path() const;
    double get_longest_dist_to(Node* dest) const;

    // Take ownership of a node
    void take_ownership(Node* node);
    // Transfer all nodes from a raw vector to owned
    void own_all(const std::vector<Node*>& nodes);
};

// Free functions for Newick parsing
std::string read_newick_from_file(const std::filesystem::path& input_file_path);
std::pair<Node*, std::vector<Node*>> root_from_newick_str(const std::string& newick_str);
std::pair<Node*, std::vector<Node*>> create_a_tree_from_newick(const std::string& newick);
Node* create_node_from_children(std::vector<std::vector<Node*>>& open_nodes_per_level,
                                 int level, double branch_length, int node_inx,
                                 std::vector<Node*>& storage);

constexpr double TREE_EPSILON = 0.001;
