#include "dist_labels_stats.h"
#include "distance_calc.h"

DistanceLabelsStats::DistanceLabelsStats(const std::string& code, int taxa_num, int msa_length)
    : BasicStats(code, taxa_num, msa_length,
                  {"code", "ssp_from_true", "dseq_from_true", "dpos_from_true"}),
      ssp_from_true(-1), dseq_from_true(-1), dpos_from_true(-1), rf_from_true(-1) {}

void DistanceLabelsStats::set_my_distance_from_true(const std::vector<std::string>& inferred_msa,
                                                    const std::vector<std::vector<std::set<std::string>>>& profile_b_h,
                                                    const std::vector<std::string>& true_msa){
    //dpos_from_true = compute_distance(true_msa, inferred_msa, DistanceType::D_POS);

    dseq_from_true = compute_distance_from_known_msa(inferred_msa, profile_b_h, DistanceType::D_SEQ);
    //ssp_from_true = compute_distance(true_msa, inferred_msa, DistanceType::D_SSP);
    dpos_from_true = compute_eff_d_seq(inferred_msa, true_msa);
    ssp_from_true = 0.0;

}

void DistanceLabelsStats::set_rf_from_true(const UnrootedTree& my_tree, const UnrootedTree& true_tree) {
    rf_from_true = my_tree.calc_rf(true_tree);
}

std::pair<std::vector<StatValue>, std::vector<std::string>> DistanceLabelsStats::get_print_rf() const {
    return {{StatValue(code), StatValue(rf_from_true)}, {"code", "rf_from_true"}};
}

std::vector<StatValue> DistanceLabelsStats::get_my_features_as_list() const {
    return {
        StatValue(code),
        StatValue(ssp_from_true), StatValue(dseq_from_true), StatValue(dpos_from_true)
    };
}
