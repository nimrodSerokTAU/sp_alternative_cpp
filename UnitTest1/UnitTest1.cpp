#include "pch.h"
#include "CppUnitTest.h"
#include "../sp_alternative_cpp/include/evo_model.h"
#include <vector>
#include "../sp_alternative_cpp/include/config.h"
#include "../sp_alternative_cpp/include/sp_score.h"
#include <numeric>
#include "../sp_alternative_cpp/include/distance_calc.h"
using namespace Microsoft::VisualStudio::CppUnitTestFramework;
using namespace std;

namespace spalternativeUnitTests
{
	TEST_CLASS(spalternativeUnitTests)
	{
	public:

		TEST_METHOD(sp_perfect)
		{
			EvoModel evoModel1(-1, -1, "Blosum50");
			vector<EvoModel> models = { evoModel1 };
			Configuration configuration1(models);
			SPScore sp1(evoModel1, "D:/code/sp_alternative/sp_alternative/input_config_files");
			vector<string> profile1 = {
				"ARNDCQEGHI",
				"ARNDCQEGHI",
				"ARNDCQEGHI"
			};
			vector<double> vw1 = { 1.0, 1.0, 1.0 };
			vector<vector<double>> vvw1 = { vw1 };
			vector<double> res = sp1.compute_naive_sp_score(profile1, &vvw1);
			// (5 + 7 + 7 + 8 + 13 + 7 + 6 + 8 + 10 + 5) * 3 = 76 * 3
			double sum = std::accumulate(res.begin(), res.end(), 0.0);
			Assert::AreEqual(sum, 228.0);
		}

		TEST_METHOD(sp_no_gaps)
		{
			EvoModel evoModel1(-1, -5, "Blosum50");
			vector<EvoModel> models = { evoModel1 };
			Configuration configuration1(models);
			SPScore sp1(evoModel1, "D:/code/sp_alternative/sp_alternative/input_config_files");
			vector<string> profile1 = {
				"ARNDCQEGHI",
				"AANDCQEGAI",
				"AANDCQEGHI"
			};
			vector<double> vw1 = { 1.0, 1.0, 1.0 };
			vector<vector<double>> vvw1 = { vw1 };
			vector<double> res = sp1.compute_naive_sp_score(profile1, &vvw1);
			//RRR->RAA : 21 -> 5 - 4 = 1 -> - 20
			//HHH->HHA : 30 -> 10 - 4 = 6 -> - 24
			double sum = std::accumulate(res.begin(), res.end(), 0.0);
			Assert::AreEqual(sum, 184.0);
		}

		TEST_METHOD(sp_local_gaps)
		{
			EvoModel evoModel1(-1, -5, "Blosum50");
			vector<EvoModel> models = { evoModel1 };
			Configuration configuration1(models);
			SPScore sp1(evoModel1, "D:/code/sp_alternative/sp_alternative/input_config_files");
			vector<string> profile1 = {
				"ARNDCQ-GHI",
				"AANDCQ-GAI",
				"AANDCQEGHI"
			};
			vector<double> vw1 = { 1.0, 1.0, 1.0 };
			vector<vector<double>> vvw1 = { vw1 };
			vector<double> res = sp1.compute_naive_sp_score(profile1, &vvw1);
			// EEE -> --E : 18 -> -6 * 2 = -12 -> -30
			double sum = std::accumulate(res.begin(), res.end(), 0.0);
			Assert::AreEqual(sum, 154.0);
		}

		TEST_METHOD(naive_algo_case_a_subs_only)
		{
			EvoModel evoModel1(0, 0, "Blosum50");
			vector<EvoModel> models = { evoModel1 };
			Configuration configuration1(models);
			SPScore sp1(evoModel1, "D:/code/sp_alternative/sp_alternative/input_config_files");
			vector<string> profile1 = {
				"ARNDC---HI",
				"AA-DCQ--AI",
				"AA--CQEGHI"
			};
			vector<double> vw1 = { 1.0, 1.0, 1.0 };
			vector<vector<double>> vvw1 = { vw1 };
			vector<double> res = sp1.compute_naive_sp_score(profile1, &vvw1);
			// 15 + 1 + (-6 - 6) + (8 - 5 - 6) + 39 + (7 - 6 - 6) + (-6 - 5) + (-5, -5) + 6 + 15
			double sum = std::accumulate(res.begin(), res.end(), 0.0);
			Assert::AreEqual(sum, 91.0);
		}

		TEST_METHOD(naive_algo_case_a_subs_and_ge)
		{
			EvoModel evoModel1(0, -5, "Blosum50");
			vector<EvoModel> models = { evoModel1 };
			Configuration configuration1(models);
			SPScore sp1(evoModel1, "D:/code/sp_alternative/sp_alternative/input_config_files");
			vector<string> profile1 = {
				"ARNDC---HI",
				"AA-DCQ--AI",
				"AA--CQEGHI"
			};
			vector<double> vw1 = { 1.0, 1.0, 1.0 };
			vector<vector<double>> vvw1 = { vw1 };
			vector<double> res = sp1.compute_naive_sp_score(profile1, &vvw1);
			// 15 + 1 + (-6 - 6) + (8 - 5 - 6) + 39 + (7 - 6 - 6) + (-6 - 5) + (-5, -5) + 6 + 15
			double sum = std::accumulate(res.begin(), res.end(), 0.0);
			Assert::AreEqual(sum, 41.0);  // ge cost only : -50
		}

		TEST_METHOD(naive_algo_case_a_subs_and_ge_and_gs)
		{
			EvoModel evoModel1(-1, -5, "Blosum50");
			vector<EvoModel> models = { evoModel1 };
			Configuration configuration1(models);
			SPScore sp1(evoModel1, "D:/code/sp_alternative/sp_alternative/input_config_files");
			vector<string> profile1 = {
				"ARNDC---HI",
				"AA-DCQ--AI",
				"AA--CQEGHI"
			};
			vector<double> vw1 = { 1.0, 1.0, 1.0 };
			vector<vector<double>> vvw1 = { vw1 };
			vector<double> res = sp1.compute_naive_sp_score(profile1, &vvw1);
			// 15 + 1 + (-6 - 6) + (8 - 5 - 6) + 39 + (7 - 6 - 6) + (-6 - 5) + (-5, -5) + 6 + 15
			double sum = std::accumulate(res.begin(), res.end(), 0.0);
			Assert::AreEqual(sum, 35.0); // gs cost only: -6
		}

		TEST_METHOD(naive_algo_case_a_subs_and_ge_and_gs_with_weights)
		{
			EvoModel evoModel1(-1, -5, "Blosum50");
			vector<EvoModel> models = { evoModel1 };
			Configuration configuration1(models);
			SPScore sp1(evoModel1, "D:/code/sp_alternative/sp_alternative/input_config_files");
			vector<string> profile1 = {
				"ARNDC---HI",
				"AA-DCQ--AI",
				"AA--CQEGHI"
			};
			vector<double> vw1 = { 1.0, 1.0, 1.0 };
			vector<double> vw2 = { 2.0, 2.0, 2.0 };
			vector<vector<double>> vvw1 = { vw1, vw2 };
			vector<double> res = sp1.compute_naive_sp_score(profile1, &vvw1);
			
			vector<double> expected = { 35.0, 140.0 };
			Assert::AreEqual(res.size(), expected.size());

			for (size_t i = 0; i < expected.size(); ++i) {
				Assert::AreEqual(res[i], expected[i]);
			}
		}

		TEST_METHOD(compute_sp_s_and_sp_ge)
		{
			EvoModel evoModel1(0, -5, "Blosum50");
			vector<EvoModel> models = { evoModel1 };
			Configuration configuration1(models);
			SPScore sp1(evoModel1, "D:/code/sp_alternative/sp_alternative/input_config_files");
			vector<string> profile1 = {
				"ARNDC---HI",
				"AA-DCQ--AI",
				"AA--CQEGHI"
			};
			SPScore::SpSAndGe res = sp1.compute_sp_s_and_sp_ge(profile1);

			vector<double> vw1 = { 1.0, 1.0, 1.0 };
			vector<vector<double>> vvw1 = { vw1 };
			vector<double> resNaive = sp1.compute_naive_sp_score(profile1, &vvw1);
			double totalScore = res.sp_match_score + res.sp_mismatch_score + res.ge_score;
			Assert::AreEqual(totalScore, resNaive[0]);
			Assert::AreEqual(res.sp_mismatch_score + res.sp_match_score, 91.0);
			Assert::AreEqual(res.ge_score, -50.0);
		}

		TEST_METHOD(only_gap_open_and_ext_cost_same)
		{
			EvoModel evoModel1(-1, -5, "Blosum50");
			vector<EvoModel> models = { evoModel1 };
			Configuration configuration1(models);
			SPScore sp1(evoModel1, "D:/code/sp_alternative/sp_alternative/input_config_files");
			vector<string> profile1 = {
				"ARNDC---HI",
				"AA-DCQ--AI",
				"AA--CQEGHI"
			};
			// ARNDC-- - HI
			// AA - DCQ--AI
			// ----> 2
			//
			// ARNDC-- - HI
			// AA--CQEGHI
			// ----> 2
			//
			// AA - DCQ--AI
			// AA--CQEGHI
			// ----> 2

			SPScore::SpGapOpen res = sp1.compute_sp_gap_open(profile1);
			Assert::AreEqual(res.sp_gp_open, -6.0);
			Assert::AreEqual(res.sp_gpo_count, 6);
		}

		TEST_METHOD(compute_efficient_sp)
		{
			EvoModel evoModel1(-1, -5, "Blosum50");
			vector<EvoModel> models = { evoModel1 };
			Configuration configuration1(models);
			SPScore sp1(evoModel1, "D:/code/sp_alternative/sp_alternative/input_config_files");
			vector<string> profile1 = {
				"ARNDC---HI",
				"AA-DCQ--AI",
				"AA--CQEGHI"
			};

			double efficientRes = sp1.compute_efficient_sp(profile1);
			vector<double> vw1 = { 1.0, 1.0, 1.0 };
			vector<vector<double>> vvw1 = { vw1 };
			vector<double> resNaive = sp1.compute_naive_sp_score(profile1, &vvw1);

			Assert::AreEqual(efficientRes, resNaive[0]);
		}

		TEST_METHOD(compare_naive_sop_to_efficient)
		{
			EvoModel evoModel1(-10, -0.5, "Blosum62");
			vector<EvoModel> models = { evoModel1 };
			Configuration configuration1(models);
			SPScore sp1(evoModel1, "D:/code/sp_alternative/sp_alternative/input_config_files");
			vector<string> profile1 = {
				"-EETTEESLKRIVADNENRAEQVHLYLSTTFVIADPEPKYGIVRSKDMNWYEQKTHKFLGMGPVLGVQFAF",
				"YEETSEESL-RIAADNENRAE-VHLYLGTNFVIADPEPKW--LRSKDVNWYDQRTH-FLGMGPVLGIQFLI",
				"YEETSEES----VADNENRAE-VHLILSTNFVIADPEPKWG-LRSKDMNWYDQRTH--LGMGPVLGIQFLF",
				"YEETSEESLKRIVADNENRAEKVHLILSTNFVIADPEPKWG--RSKDMNWYDQRTHKFLGMGPVLGIQFLF"
			};

			double efficientRes = sp1.compute_efficient_sp(profile1);
			vector<double> vw1 = { 1.0, 1.0, 1.0, 1.0 };
			vector<vector<double>> vvw1 = { vw1 };
			vector<double> resNaive = sp1.compute_naive_sp_score(profile1, &vvw1);

			Assert::AreEqual(efficientRes, resNaive[0]);
		}

		TEST_METHOD(compare_naive_sop_to_efficient_min_example)
		{
			EvoModel evoModel1(-10, -0.5, "Blosum62");
			vector<EvoModel> models = { evoModel1 };
			Configuration configuration1(models);
			SPScore sp1(evoModel1, "D:/code/sp_alternative/sp_alternative/input_config_files");
			vector<string> profile1 = {
				"LLKYR-K",
				"Y--ERAK",
				"YL----K",
				"YLKE-AK"
			};

			double efficientRes = sp1.compute_efficient_sp(profile1);
			vector<double> vw1 = { 1.0, 1.0, 1.0, 1.0 };
			vector<vector<double>> vvw1 = { vw1 };
			vector<double> resNaive = sp1.compute_naive_sp_score(profile1, &vvw1);

			Assert::AreEqual(efficientRes, resNaive[0]);
		}

		TEST_METHOD(translate_profile_hpos)
		{
			vector<string> profile1 = {
				"AATATTG-",
				"A--ATTAG",
				"A--A-TAG"
			};

			vector<vector<string>> res = translate_profile_naming(profile1, DistanceType::D_POS);
			vector<string> vw1 = { "S^1_1", "S^1_2", "S^1_3", "S^1_4", "S^1_5", "S^1_6", "S^1_7", "G^1_7" };
			vector<string> vw2 = { "S^2_1", "G^2_1", "G^2_1", "S^2_2", "S^2_3", "S^2_4", "S^2_5", "S^2_6" };
			vector<string> vw3 = { "S^3_1", "G^3_1", "G^3_1", "S^3_2", "G^3_2", "S^3_3", "S^3_4", "S^3_5" };
			vector<vector<string>> expected = { vw1, vw2, vw3 };

			Assert::AreEqual(expected.size(), res.size());

			for (size_t i = 0; i < expected.size(); ++i)
			{
				Assert::AreEqual(expected[i].size(), res[i].size());

				for (size_t j = 0; j < expected[i].size(); ++j)
				{
					Assert::AreEqual(
						std::wstring(expected[i][j].begin(), expected[i][j].end()),
						std::wstring(res[i][j].begin(), res[i][j].end())
					);
				}
			}

		}

		TEST_METHOD(get_specific_hpos)
		{ 
			vector<string> profile1 = {
				"AATATTG-",
				"A--ATTAG",
				"A--A-TAG"
			};

			vector<vector<string>> trans_prof = translate_profile_naming(profile1, DistanceType::D_POS);

			std::vector<std::optional<std::set<std::string>>> res_2;
			std::vector<std::optional<std::set<std::string>>> res_3;

			for (size_t col_inx = 0; col_inx < profile1[0].size(); ++col_inx)
			{
				std::vector<std::string> col = get_column(trans_prof, col_inx);

				std::optional<std::set<std::string>> hpos_2 = get_place_h(col, 1);
				res_2.push_back(hpos_2);

				std::optional<std::set<std::string>> hpos_3 = get_place_h(col, 2);
				res_3.push_back(hpos_3);
			}

			std::optional<std::set<std::string>> res_2_0 = std::set<std::string>{ "S^1_1", "S^3_1" };
			std::optional<std::set<std::string>> res_2_3 = std::set<std::string>{ "S^1_4", "S^3_2" };
			std::optional<std::set<std::string>> res_2_4 = std::set<std::string>{ "S^1_5", "G^3_2" };
			std::optional<std::set<std::string>> res_2_5 = std::set<std::string>{ "S^1_6", "S^3_3" };
			std::optional<std::set<std::string>> res_2_6 = std::set<std::string>{ "S^1_7", "S^3_4" };
			std::optional<std::set<std::string>> res_2_7 = std::set<std::string>{ "G^1_7", "S^3_5" };

			vector<std::optional<std::set<std::string>>> res2_exprected = { res_2_0, std::nullopt, std::nullopt, res_2_3, res_2_4, res_2_5, res_2_6, res_2_7 };

			Assert::AreEqual(res2_exprected.size(), res_2.size());

			for (size_t i = 0; i < res2_exprected.size(); ++i)
			{
				Assert::IsTrue(res2_exprected[i] == res_2[i]);
			}

			std::optional<std::set<std::string>> res_3_0 = std::set<std::string>{ "S^1_1", "S^2_1" };
			std::optional<std::set<std::string>> res_3_3 = std::set<std::string>{ "S^1_4", "S^2_2" };
			std::optional<std::set<std::string>> res_3_5 = std::set<std::string>{ "S^1_6", "S^2_4" };
			std::optional<std::set<std::string>> res_3_6 = std::set<std::string>{ "S^1_7", "S^2_5" };
			std::optional<std::set<std::string>> res_3_7 = std::set<std::string>{ "G^1_7", "S^2_6" };

			vector<std::optional<std::set<std::string>>> res3_exprected = { res_3_0, std::nullopt, std::nullopt, res_3_3, std::nullopt, res_3_5, res_3_6, res_3_7 };

			Assert::AreEqual(res3_exprected.size(), res_3.size());

			for (size_t i = 0; i < res3_exprected.size(); ++i)
			{
				Assert::IsTrue(res3_exprected[i] == res_3[i]);
			}

		}

		TEST_METHOD(compute_dpos_distance_for_same)
		{
			vector<string> profile1 = {
				"AATATTG-",
				"A--ATTAG",
				"A--A-TAG"
			};
			double res = compute_distance(profile1, profile1, DistanceType::D_SEQ);
			Assert::AreEqual(res, 0.0);
		}

		TEST_METHOD(compute_dpos_distance_for_diff)
		{
			vector<string> profile1 = {
				"AATATTG-",
				"A--ATTAG",
				"A--A-TAG"
			};

			vector<string> profile2 = {
				"AATAT-TG",
				"A-A-TTAG",
				"A--A-TAG"
			};
			double res = compute_distance(profile1, profile2, DistanceType::D_POS);
			Assert::AreEqual(res, 0.417, 0.001);
		}

		TEST_METHOD(dpos_distance_for_diff_case_a)
		{
			vector<string> profile1 = {
				"AATATTG-",
				"A--ATTAG",
				"A--A-TAG"
			};

			vector<string> profile2 = {
				"A-ATAT-TG",
				"A-A--TTAG",
				"A--A-T-AG"
			};
			double res = compute_distance(profile1, profile2, DistanceType::D_POS);
			Assert::AreEqual(res, 0.639, 0.001);
		}

		TEST_METHOD(dpos_for_diff_length_case_qu_a)
		{
			vector<string> profile1 = {
				"GCATCATT-G",
				"GC---ATTAG",
				"GC---AT-AG"
			};

			vector<string> profile2 = {
				"GCATCATT-G",
				"GCA---TTAG",
				"GCA----TAG"
			};

			vector<string> profile3 = {
				"GCATCATT-G-",
				"GC---ATT-AG",
				"GC---A-TA-G"
			};

			double res_1_2 = compute_distance(profile1, profile2, DistanceType::D_POS);
			double res_1_3 = compute_distance(profile1, profile3, DistanceType::D_POS);
			Assert::AreEqual(res_1_2, 0.364, 0.001);
			Assert::AreEqual(res_1_3, 0.295, 0.001);
		}

		TEST_METHOD(dseq_for_diff_length_case_qu_a)
		{
			vector<string> profile1 = {
				"GCATCATT-G",
				"GC---ATTAG",
				"GC---AT-AG"
			};

			vector<string> profile2 = {
				"GCATCATT-G",
				"GCA---TTAG",
				"GCA----TAG"
			};

			vector<string> profile3 = {
				"GCATCATT-G-",
				"GC---ATT-AG",
				"GC---A-TA-G"
			};

			double res_1_2 = compute_distance(profile1, profile2, DistanceType::D_SEQ);
			double res_1_3 = compute_distance(profile1, profile3, DistanceType::D_SEQ);
			Assert::AreEqual(res_1_2, 0.273, 0.001);
			Assert::AreEqual(res_1_3, 0.265, 0.001);
		}

		TEST_METHOD(dssp_for_diff_length_case_qu_a)
		{
			vector<string> profile1 = {
				"GCATCATT-G",
				"GC---ATTAG",
				"GC---AT-AG"
			};

			vector<string> profile2 = {
				"GCATCATT-G",
				"GCA---TTAG",
				"GCA----TAG"
			};

			vector<string> profile3 = {
				"GCATCATT-G-",
				"GC---ATT-AG",
				"GC---A-TA-G"
			};

			double res_1_2 = compute_distance(profile1, profile2, DistanceType::D_SSP);
			double res_1_3 = compute_distance(profile1, profile3, DistanceType::D_SSP);
			Assert::AreEqual(res_1_2, 0.381, 0.001);
			Assert::AreEqual(res_1_3, 0.4, 0.001);
		}


	};
}
