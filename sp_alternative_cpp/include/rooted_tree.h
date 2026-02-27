#pragma once
#include "node.h"
#include "unrooted_tree.h"
#include "enums.h"
#include <vector>
#include <unordered_map>
#include <string>
#include <deque>

struct RootingPoint {
    int start_id;
    int end_id;
    double dist_from_start;
    double dist_from_end;
    double longest_dist; // used for shallowest tree comparison
};

class RootedTree {
public:
    Node* root;
    std::vector<Node*> all_nodes;
    std::set<std::string> keys;
    std::unordered_map<std::string, double> seq_weight_dict;

    // Memory management for deep-copied nodes
    std::vector<std::unique_ptr<Node>> owned_nodes;

    RootedTree(Node* root, const std::vector<Node*>& all_nodes, const std::set<std::string>& keys);

    static RootedTree root_tree(const UnrootedTree& unrooted, RootingMethods rooting_method);

    void calc_clustal_w();
    void calc_seq_w();
};

// Helper free functions
RootingPoint calc_mid_point(const UnrootedTree& unrooted);
std::vector<RootingPoint> calc_min_differential_sum(const UnrootedTree& unrooted,
                                                      std::vector<Node*>& all_nodes);
RootingPoint find_shallowest_tree(std::vector<Node*>& all_nodes,
                                   std::vector<RootingPoint>& rooting_points);
void recalc_tree_down(Node* node, Node* father, int broke_id,
                       std::vector<std::tuple<Node*, Node*, int>>& nodes_to_recalc,
                       std::vector<Node*>& all_nodes);
void fill_nodes_w(Node* node, std::deque<Node*>& nodes_to_recalc);
double sum_bl_up_to_node_id(Node* origin, int dest_id, const std::vector<Node*>& all_nodes);
RootingPoint calc_potential_root_on_branch(double bl, double w_to_orig, double w_to_dest,
                                            double delta, int origin, int dest, double min_bl);
std::pair<Node*, std::vector<Node*>> create_root(std::vector<Node*>& all_nodes,
                                                   const RootingPoint& rooting_point,
                                                   std::vector<std::unique_ptr<Node>>& storage);

// Deep copy utilities
std::vector<Node*> deep_copy_nodes(const std::vector<Node*>& src,
                                    std::vector<std::unique_ptr<Node>>& storage);
