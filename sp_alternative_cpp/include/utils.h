#pragma once
#include <vector>
#include <string>
#include <unordered_map>
#include <utility>

std::pair<std::vector<std::vector<int>>, std::unordered_map<std::string, int>>
read_matching_matrix(const std::string& file_path);

int translate_to_matrix_index(char letter, const std::unordered_map<std::string, int>& code_to_index_dict);

double calc_p_distance_from_other(const std::string& aligned_seq, const std::string& other_aligned_seq);

double calc_kimura_distance_from_other(const std::string& aligned_seq, const std::string& other_aligned_seq);

double calc_percentile(std::vector<double>& values, int percentile);
double calc_percentile_int(std::vector<int>& values, int percentile);
