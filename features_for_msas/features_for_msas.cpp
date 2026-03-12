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

	set<string> weightMethodsSet = config["additional_weights"].get<set<string>>();
    set<WeightMethods> resultWeightMethodsSet = strToEnumWeightMethods(weightMethodsSet);

    set<string> statsOutputSet = config["stats_output"].get<set<string>>();
    set<StatsOutput> resultStatsOutputSet = strToEnumStatsOutput(statsOutputSet);

    set<int> k_values = config["k_values"].get<set<int>>();

    int x = config["sop_calc_type"].get<int>();

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