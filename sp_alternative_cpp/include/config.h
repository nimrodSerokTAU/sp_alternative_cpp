#pragma once
#include <vector>
#include <set>
#include <string>
#include "evo_model.h"
#include "enums.h"

struct Configuration {
    std::vector<EvoModel> models;
    SopCalcTypes sop_calc_type;
    std::string input_files_dir_name;
    std::set<WeightMethods> additional_weights;
    std::set<StatsOutput> stats_output;
    std::set<int> k_values;

    Configuration(const std::vector<EvoModel>& models_list,
                  SopCalcTypes sop_calc_type = SopCalcTypes::EFFICIENT,
                  const std::string& input_files_dir_name = "",
                  const std::set<WeightMethods>& additional_weights = {},
                  const std::set<int>& k_values = {},
                  const std::set<StatsOutput>& stats_output = {StatsOutput::ALL})
        : models(models_list),
          sop_calc_type(sop_calc_type),
          input_files_dir_name(input_files_dir_name),
          additional_weights(additional_weights),
          stats_output(stats_output),
          k_values(k_values) {}
};
