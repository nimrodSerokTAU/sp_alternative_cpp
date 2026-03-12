#include <iostream>
#include <fstream>
#include <nlohmann/json.hpp>
#include <string>
#include <vector>
#include <tuple>
#include <optional>
#include <filesystem>
#include <unordered_map>
#include <chrono>
using namespace std;

#include "../sp_alternative_cpp/include/evo_model.h"
#include "../sp_alternative_cpp/include/config.h"
#include "../sp_alternative_cpp/include/sp_score.h"
#include "../sp_alternative_cpp/include/msa.h"
#include "../sp_alternative_cpp/include/multi_msa_service.h"
#include "../sp_alternative_cpp/include/enums.h"

using json = nlohmann::json;

int main()
{
    std::ifstream f("../config.json");
    json config = json::parse(f);
    std::vector<EvoModel> models;
    for (const auto& modelConfig : config["models_list"]) {
        int gapOpenCost = modelConfig["gap_open_cost"].get<int>();
        int gapExtendCost = modelConfig["gap_extend_cost"].get<int>();
        std::string substitutionMatrix = modelConfig["matrix_file_name"].get<std::string>();
        EvoModel evoModel(gapOpenCost, gapExtendCost, substitutionMatrix);
        models.push_back(evoModel);
    }

	int x = config["sop_calc_type"].get<int>();
	
    set<string> weightMethodsSet = config["additional_weights"].get<set<string>>();

    set<WeightMethods> resultWeightMethodsSet;

    for (const auto& s : weightMethodsSet) {
        auto it = str_to_enum_weight_methods.find(s);
        if (it != str_to_enum_weight_methods.end()) {
            resultWeightMethodsSet.insert(it->second);
        }
    }

    set<string> statsOutputSet = config["stats_output"].get<set<string>>();

    set<StatsOutput> resultStatsOutputSet;

    for (const auto& s : statsOutputSet) {
        auto it = str_to_enum_stats_output.find(s);
        if (it != str_to_enum_stats_output.end()) {
            resultStatsOutputSet.insert(it->second);
        }
    }

    set<int> k_values = config["k_values"].get<set<int>>();

    Configuration configuration(
        models,
        static_cast<SopCalcTypes>(x),
        config["input_files_dir_path"],
        config["output_file_dir_path"],
        config["matrix_dir_path"],
        resultWeightMethodsSet,
        k_values,
        resultStatsOutputSet
	);
    cout << "SP Alternative (C++) - MSA Analysis Tool" << endl;
    cout << "Input directory: " << configuration.input_files_dir_path << endl;

    try {
        multiple_msa_calc_features_and_labels(configuration);
    }
    catch (const exception& e) {
        cerr << "Error: " << e.what() << endl;
        return 1;
    }

    return 0;;
}

//std::tuple<std::optional<std::string>, std::optional<std::string>, std::vector<std::string>>
//get_file_names_ordered(const std::vector<std::string>& file_names)
//{
//    std::optional<std::string> true_file_name;
//    std::optional<std::string> true_tree_file_name;
//    std::vector<std::string> other_file_names;
//
//    for (const auto& file_name : file_names) {
//        if (file_name.find("_TRUE.") != std::string::npos) {
//            true_file_name = file_name;
//        } else if (file_name.find("_true_tree.") != std::string::npos) {
//            true_tree_file_name = file_name;
//        } else {
//            other_file_names.push_back(file_name);
//        }
//    }
//
//    return { true_file_name, true_tree_file_name, other_file_names };
//}
//
//namespace fs = std::filesystem;
//
//void multiple_msa_calc_features_and_labels(const Configuration& config)
//{
//    auto start = std::chrono::steady_clock::now();
//
//    std::unordered_map<std::string, std::unordered_map<std::string, double>> all_msa_stats;
//
//    for (const auto& stats_file_name : config.stats_output) {
//        all_msa_stats[stats_output_to_string(stats_file_name)] = {};
//    }
//
//    fs::path comparison_dir = config.input_files_dir_path;
//
//    std::vector<SPScore> sp_models;
//    for (const auto& m : config.models) {
//        sp_models.emplace_back(m, "");
//    }
//
//    fs::path output_dir_path =  config.output_file_dir_path;
//
//    for (const auto& dir_entry : fs::directory_iterator(comparison_dir)) {
//        if (!dir_entry.is_directory())
//            continue;
//
//        fs::path dir_path = dir_entry.path();
//        std::string dir_name = dir_path.filename().string();
//
//        std::vector<std::string> file_names;
//        for (const auto& file : fs::directory_iterator(dir_path)) {
//            file_names.push_back(file.path().filename().string());
//        }
//
//        auto [true_file_name, true_tree_file_name, inferred_file_names] =
//            get_file_names_ordered(file_names);
//
//        MSA true_msa(dir_name);
//
//        if (true_file_name.has_value()) {
//            true_msa.read_from_fasta(dir_path / true_file_name.value());
//        }
//
//        bool need_tree =
//            config.stats_output.count(StatsOutput::ALL) ||
//            config.stats_output.count(StatsOutput::RF_LABEL);
//
//        if (need_tree) {
//            if (true_tree_file_name.has_value()) {
//                true_msa.set_tree(
//                    UnrootedTree::create_from_newick_file(
//                        dir_path / true_tree_file_name.value()
//                    )
//                );
//            } else {
//                true_msa.build_nj_tree();
//            }
//        }
//
//        bool is_init_files = true;
//
//        for (const auto& inferred_file_name : inferred_file_names) {
//            std::string msa_name = inferred_file_name;
//            std::cout << msa_name << std::endl;
//
//            MSA inferred_msa(msa_name);
//            inferred_msa.read_from_fasta(dir_path / inferred_file_name);
//            inferred_msa.order_sequences(true_msa.seq_names);
//
//            inferred_msa.calc_and_print_stats(
//                true_msa,
//                config,
//                sp_models,
//                output_dir_path,
//                true_msa.tree ? true_msa.tree.get() : nullptr,
//                is_init_files
//            );
//
//            is_init_files = false;
//        }
//
//        true_msa.calc_and_print_stats(
//            true_msa,
//            config,
//            sp_models,
//            output_dir_path,
//            true_msa.tree ? true_msa.tree.get() : nullptr,
//            is_init_files
//        );
//    }
//
//    auto end = std::chrono::steady_clock::now();
//    auto duration =
//        std::chrono::duration_cast<std::chrono::seconds>(end - start).count();
//
//    std::cout << "Done - Total duration: " << duration << " seconds" << std::endl;
//}
