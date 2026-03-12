#include "pch.h"
#include "CppUnitTest.h"
#include "../sp_alternative_cpp/include/evo_model.h"
#include <vector>
#include "../sp_alternative_cpp/include/config.h"
#include "../sp_alternative_cpp/include/sp_score.h"
#include <numeric>
#include "../sp_alternative_cpp/include/distance_calc.h"
#include "../sp_alternative_cpp/include/node.h"
#include "../sp_alternative_cpp/include/unrooted_tree.h"
using namespace Microsoft::VisualStudio::CppUnitTestFramework;
using namespace std;

namespace spalternativeUnitTests
{
	string blosum50Path = "D:/code/sp_alternative/sp_alternative/input_config_files/Blosum50.txt";
	string blosum62Path = "D:/code/sp_alternative/sp_alternative/input_config_files/Blosum62.txt";
	string newick_of_AATF = "((((Macropus:0.051803,Monodelphis:  0.066021): 0.016682,Sarcophilus:0.068964):0.114355,((Echinops:0.104144,(Loxodonta:0.076474,Procavia:0.076193):0.011550):0.013015,((Choloepus:0.056091,Dasypus:0.040600):0.013681,(((((Callithrix:0.032131,((((Gorilla:0.007042,(Homo:0.002445,Pan:0.002450):0.001237):0.003689,Pongo:0.007508):0.002384,Nomascus:0.015696):0.004674,Macaca:0.013752):0.012205):0.029730,((Microcebus:0.037066,Otolemur:0.050935):0.008091,Tarsius:0.064938):0.007305):0.003377,Tupaia:0.090946):0.000699,(((Cavia:0.116826,(Dipodomys:0.080386,(Mus:0.040313,Rattus:0.033872):0.122329):0.013314):0.000298,Ictidomys:0.062932):0.011064,(Ochotona:0.087746,Oryctolagus:0.057769):0.035947):0.002130):0.006931,(Erinaceus:0.094727,(((((Bos:0.062659,Tursiops:0.024374):0.009782,Sus:0.064336):0.006134,Vicugna:0.049961):0.021643,Equus:0.046559):0.002728,(Sorex:0.126677,((Myotis:0.050452,Pteropus:0.047740):0.006495,(Felis:0.042414,(Canis:0.036675,(Mustela:0.027691,Ailuropoda:0.037553):0.004251):0.008141):0.019485):0.000485):0.003533):0.000528):0.005645):0.012900):0.002853):0.133251):0.019558,Ornithorhynchus:0.195576);";

	TEST_CLASS(spalternativeUnitTests)
	{
	public:

		TEST_METHOD(sp_perfect)
		{
			EvoModel evoModel1(-1, -1, "Blosum50");
			vector<EvoModel> models = { evoModel1 };
			Configuration configuration1(models);
			SPScore sp1(evoModel1, blosum50Path);
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
			SPScore sp1(evoModel1, blosum50Path);
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
			SPScore sp1(evoModel1, blosum50Path);
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
			SPScore sp1(evoModel1, blosum50Path);
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
			SPScore sp1(evoModel1, blosum50Path);
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
			SPScore sp1(evoModel1, blosum50Path);
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
			SPScore sp1(evoModel1, blosum50Path);
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
			SPScore sp1(evoModel1, blosum50Path);
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
			SPScore sp1(evoModel1, blosum50Path);
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
			SPScore sp1(evoModel1, blosum50Path);
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
			SPScore sp1(evoModel1, blosum62Path);
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
			SPScore sp1(evoModel1, blosum62Path);
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

			for (int col_inx = 0; col_inx < profile1[0].size(); ++col_inx)
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

		TEST_METHOD(dpos_for_diff_length_case_qu_b)
		{
			vector<string> trueProfile = {
				"GAAGTTAGACATC",
				"GA--------ATC",
				"GA--------ATG",
				"GT--------ATG",
				"GT--------ATC"
			};

			vector<string> profileA = {
				"GAAGTTAGACATC",
				"GAA--------TC",
				"GAA--------TG",
				"GTA--------TG",
				"GTA--------TC"
			};

			vector<string> profileB = {
				"GAAGTTAGACATC",
				"GA----A----TC",
				"GA----A----TG",
				"GT----A----TG",
				"GT----A----TC"
			};

			vector<string> profileC = {
				"GAAGTTAGACATC",
				"GA------A--TC",
				"GA------A--TG",
				"GT------A--TG",
				"GT------A--TC"
			};

			double res_t_1 = compute_distance(trueProfile, profileA, DistanceType::D_POS);
			double res_t_2 = compute_distance(trueProfile, profileB, DistanceType::D_POS);
			double res_t_3 = compute_distance(trueProfile, profileC, DistanceType::D_POS);
			Assert::AreEqual(res_t_1, 0.303, 0.001);
			Assert::AreEqual(res_t_2, 0.182, 0.001);
			Assert::AreEqual(res_t_3, 0.121, 0.001);
		}

		TEST_METHOD(dseq_for_diff_length_case_qu_b)
		{
			vector<string> trueProfile = {
				"GAAGTTAGACATC",
				"GA--------ATC",
				"GA--------ATG",
				"GT--------ATG",
				"GT--------ATC"
			};

			vector<string> profileA = {
				"GAAGTTAGACATC",
				"GAA--------TC",
				"GAA--------TG",
				"GTA--------TG",
				"GTA--------TC"
			};

			vector<string> profileB = {
				"GAAGTTAGACATC",
				"GA----A----TC",
				"GA----A----TG",
				"GT----A----TG",
				"GT----A----TC"
			};

			vector<string> profileC = {
				"GAAGTTAGACATC",
				"GA------A--TC",
				"GA------A--TG",
				"GT------A--TG",
				"GT------A--TC"
			};

			double res_t_1 = compute_distance(trueProfile, profileA, DistanceType::D_SEQ);
			double res_t_2 = compute_distance(trueProfile, profileB, DistanceType::D_SEQ);
			double res_t_3 = compute_distance(trueProfile, profileC, DistanceType::D_SEQ);
			Assert::AreEqual(res_t_1, 0.091, 0.001);
			Assert::AreEqual(res_t_2, 0.091, 0.001);
			Assert::AreEqual(res_t_3, 0.091, 0.001);
		}

		TEST_METHOD(dssp_for_diff_length_case_qu_b)
		{
			vector<string> trueProfile = {
				"GAAGTTAGACATC",
				"GA--------ATC",
				"GA--------ATG",
				"GT--------ATG",
				"GT--------ATC"
			};

			vector<string> profileA = {
				"GAAGTTAGACATC",
				"GAA--------TC",
				"GAA--------TG",
				"GTA--------TG",
				"GTA--------TC"
			};

			vector<string> profileB = {
				"GAAGTTAGACATC",
				"GA----A----TC",
				"GA----A----TG",
				"GT----A----TG",
				"GT----A----TC"
			};

			vector<string> profileC = {
				"GAAGTTAGACATC",
				"GA------A--TC",
				"GA------A--TG",
				"GT------A--TG",
				"GT------A--TC"
			};

			double res_t_1 = compute_distance(trueProfile, profileA, DistanceType::D_SSP);
			double res_t_2 = compute_distance(trueProfile, profileB, DistanceType::D_SSP);
			double res_t_3 = compute_distance(trueProfile, profileC, DistanceType::D_SSP);
			Assert::AreEqual(res_t_1, 0.148, 0.001);
			Assert::AreEqual(res_t_2, 0.148, 0.001);
			Assert::AreEqual(res_t_3, 0.148, 0.001);
		}

		TEST_METHOD(dpos_for_diff_length_case_c)
		{
			vector<string> profileA = {
				"AAT",
				"--T",
				"-CG"
			};

			vector<string> profileB = {
				"AAT-",
				"---T",
				"--CG"
			};

			double res = compute_distance(profileB, profileA, DistanceType::D_POS);
			Assert::AreEqual(res, 0.5);

		}

		TEST_METHOD(dpos_for_diff_length_case_d)
		{
			vector<string> profileA = {
				"ATA",
				"-T-",
				"CG-"
			};

			vector<string> profileB = {
				"AT-A",
				"--T-",
				"-CG-"
			};

			double res = compute_distance(profileB, profileA, DistanceType::D_POS);
			Assert::AreEqual(res, 0.5);
		}

		TEST_METHOD(tree_from_newick)
		{
			std::pair<Node*, std::vector<Node*>> res = create_a_tree_from_newick(newick_of_AATF);
			std::unordered_map<int, std::set<std::string>> expextedTaxaDict;
			std::unordered_map<int, std::vector<int>> expextedChildrenIdsDict;

			expextedChildrenIdsDict[78] = {76, 77};
			expextedChildrenIdsDict[76] = {4, 75};
			expextedTaxaDict[76] = {
				"Ailuropoda", "Bos", "Callithrix", "Canis", "Cavia", "Choloepus", "Dasypus", "Dipodomys", "Echinops",
				"Equus", "Erinaceus", "Felis", "Gorilla", "Homo", "Ictidomys", "Loxodonta", "Macaca", "Macropus",
				"Microcebus", "Monodelphis", "Mus", "Mustela", "Myotis", "Nomascus", "Ochotona", "Oryctolagus", "Otolemur",	"Pan", 
				"Pongo", "Procavia", "Pteropus", "Rattus", "Sarcophilus", "Sorex", "Sus", "Tarsius", "Tupaia", "Tursiops", "Vicugna" };
			expextedChildrenIdsDict[75] = { 9, 74 };
			expextedTaxaDict[75] = {
				"Ailuropoda", "Bos", "Callithrix", "Canis", "Cavia", "Choloepus", "Dasypus", "Dipodomys", "Echinops",
				"Equus", "Erinaceus", "Felis", "Gorilla", "Homo", "Ictidomys", "Loxodonta", "Macaca", "Microcebus",
				"Mus", "Mustela", "Myotis", "Nomascus", "Ochotona", "Oryctolagus", "Otolemur", "Pan", "Pongo", "Procavia",
				"Pteropus", "Rattus", "Sorex", "Sus", "Tarsius", "Tupaia", "Tursiops", "Vicugna" };
			expextedChildrenIdsDict[74] = { 12, 73 };
			expextedTaxaDict[74] = {
				"Ailuropoda", "Bos", "Callithrix", "Canis", "Cavia", "Choloepus", "Dasypus", "Dipodomys", "Equus","Erinaceus", "Felis", "Gorilla", 
				"Homo", "Ictidomys", "Macaca", "Microcebus", "Mus", "Mustela", "Myotis", "Nomascus", "Ochotona", "Oryctolagus", "Otolemur", "Pan", 
				"Pongo", "Pteropus", "Rattus", "Sorex", "Sus", "Tarsius", "Tupaia", "Tursiops", "Vicugna" };
			expextedChildrenIdsDict[73] = { 47, 72 };
			expextedTaxaDict[73] = { "Ailuropoda", "Bos", "Callithrix", "Canis", "Cavia", "Dipodomys", "Equus", "Erinaceus", "Felis", "Gorilla",
				"Homo", "Ictidomys", "Macaca", "Microcebus", "Mus", "Mustela", "Myotis", "Nomascus", "Ochotona", "Oryctolagus", "Otolemur", "Pan", "Pongo", 
				"Pteropus", "Rattus", "Sorex", "Sus", "Tarsius", "Tupaia", "Tursiops", "Vicugna" };
			expextedChildrenIdsDict[47] = { 33, 46 };
			expextedTaxaDict[47] = { "Callithrix", "Cavia", "Dipodomys", "Gorilla", "Homo", "Ictidomys", "Macaca", "Microcebus", "Mus",
				"Nomascus", "Ochotona", "Oryctolagus", "Otolemur", "Pan", "Pongo", "Rattus", "Tarsius", "Tupaia" };
			expextedChildrenIdsDict[72] = { 48, 71 };
			expextedTaxaDict[72] = { "Ailuropoda", "Bos", "Canis", "Equus", "Erinaceus", "Felis", "Mustela", "Myotis", "Pteropus", "Sorex",
				"Sus", "Tursiops", "Vicugna" };
			expextedChildrenIdsDict[71] = { 57, 70 };
			expextedTaxaDict[71] = { "Ailuropoda", "Bos", "Canis", "Equus", "Felis", "Mustela", "Myotis", "Pteropus", "Sorex", "Sus", "Tursiops", "Vicugna" };
			expextedChildrenIdsDict[33] = { 31, 32 };
			expextedTaxaDict[33] = { "Callithrix", "Gorilla", "Homo", "Macaca", "Microcebus", "Nomascus", "Otolemur", "Pan", "Pongo", "Tarsius", "Tupaia" };
			expextedChildrenIdsDict[31] = { 25, 30 };
			expextedTaxaDict[31] = { "Callithrix", "Gorilla", "Homo", "Macaca", "Microcebus", "Nomascus", "Otolemur", "Pan", "Pongo", "Tarsius" };
			expextedChildrenIdsDict[70] = { 58, 69 };
			expextedTaxaDict[70] = { "Ailuropoda", "Canis", "Felis", "Mustela", "Myotis", "Pteropus", "Sorex" };
			expextedChildrenIdsDict[46] = { 42, 45 };
			expextedTaxaDict[46] = { "Cavia", "Dipodomys", "Ictidomys", "Mus", "Ochotona", "Oryctolagus", "Rattus" };
			expextedChildrenIdsDict[25] = { 13, 24 };
			expextedTaxaDict[25] = { "Callithrix", "Gorilla", "Homo", "Macaca", "Nomascus", "Pan", "Pongo" };
			expextedChildrenIdsDict[69] = { 61, 68 };
			expextedTaxaDict[69] = { "Pteropus", "Mustela", "Canis", "Felis", "Ailuropoda", "Myotis" };
			expextedChildrenIdsDict[24] = { 22, 23 };
			expextedTaxaDict[24] = { "Pan", "Gorilla", "Macaca", "Homo", "Pongo", "Nomascus" };
			expextedChildrenIdsDict[57] = { 55, 56 };
			expextedTaxaDict[57] = { "Vicugna", "Sus", "Equus", "Bos", "Tursiops" };
			expextedChildrenIdsDict[42] = { 40, 41 };
			expextedTaxaDict[42] = { "Cavia", "Mus", "Dipodomys", "Ictidomys", "Rattus" };
			expextedChildrenIdsDict[22] = { 20, 21 };
			expextedTaxaDict[22] = { "Pan", "Gorilla", "Homo", "Pongo", "Nomascus" };
			expextedChildrenIdsDict[68] = { 62, 67 };
			expextedTaxaDict[68] = { "Canis", "Felis", "Mustela", "Ailuropoda" };
			expextedChildrenIdsDict[55] = { 53, 54 };
			expextedTaxaDict[55] = { "Sus", "Vicugna", "Bos", "Tursiops" };
			expextedChildrenIdsDict[40] = { 34, 39 };
			expextedTaxaDict[40] = { "Dipodomys", "Cavia", "Rattus", "Mus" };
			expextedChildrenIdsDict[20] = { 18, 19 };
			expextedTaxaDict[20] = { "Homo", "Pongo", "Pan", "Gorilla" };
			expextedChildrenIdsDict[67] = { 63, 66 };
			expextedTaxaDict[67] = { "Canis", "Mustela", "Ailuropoda" };
			expextedChildrenIdsDict[53] = { 51, 52 };
			expextedTaxaDict[53] = { "Sus", "Bos", "Tursiops" };
			expextedChildrenIdsDict[39] = { 35, 38 };
			expextedTaxaDict[39] = { "Dipodomys", "Rattus", "Mus" };
			expextedChildrenIdsDict[30] = { 28, 29 };
			expextedTaxaDict[30] = { "Microcebus", "Otolemur", "Tarsius" };
			expextedChildrenIdsDict[18] = { 14, 17 };
			expextedTaxaDict[18] = { "Homo", "Pan", "Gorilla" };
			expextedChildrenIdsDict[9] = { 5 ,  8 };
			expextedTaxaDict[9] = { "Procavia", "Loxodonta", "Echinops" };
			expextedChildrenIdsDict[4] = { 2 ,  3 };
			expextedTaxaDict[4] = { "Monodelphis", "Macropus", "Sarcophilus" };
			expextedChildrenIdsDict[66] = { 64, 65 };
			expextedTaxaDict[66] = { "Mustela", "Ailuropoda" };
			expextedChildrenIdsDict[61] = { 59, 60 };
			expextedTaxaDict[61] = { "Pteropus", "Myotis" };
			expextedChildrenIdsDict[51] = { 49, 50 };
			expextedTaxaDict[51] = { "Bos", "Tursiops" };
			expextedChildrenIdsDict[45] = { 43, 44 };
			expextedTaxaDict[45] = { "Oryctolagus", "Ochotona" };
			expextedChildrenIdsDict[38] = { 36, 37 };
			expextedTaxaDict[38] = { "Rattus", "Mus" };
			expextedChildrenIdsDict[28] = { 26, 27 };
			expextedTaxaDict[28] = { "Microcebus", "Otolemur" };
			expextedChildrenIdsDict[17] = { 15, 16 };
			expextedTaxaDict[17] = { "Homo", "Pan" };
			expextedChildrenIdsDict[12] = { 10, 11 };
			expextedTaxaDict[12] = { "Dasypus", "Choloepus" };
			expextedChildrenIdsDict[8] = { 6, 7 };
			expextedTaxaDict[8] = { "Procavia", "Loxodonta" };
			expextedChildrenIdsDict[2] = { 0, 1 };
			expextedTaxaDict[2] = { "Monodelphis", "Macropus" };
			expextedTaxaDict[77] = { "Ornithorhynchus" };
			expextedTaxaDict[65] = { "Ailuropoda" };
			expextedTaxaDict[64] = { "Mustela" };
			expextedTaxaDict[63] = { "Canis" };
			expextedTaxaDict[62] = { "Felis" };
			expextedTaxaDict[60] = { "Pteropus" };
			expextedTaxaDict[59] = { "Myotis" };
			expextedTaxaDict[58] = { "Sorex" };
			expextedTaxaDict[56] = { "Equus" };
			expextedTaxaDict[54] = { "Vicugna" };
			expextedTaxaDict[52] = { "Sus" };
			expextedTaxaDict[50] = { "Tursiops" };
			expextedTaxaDict[49] = { "Bos" };
			expextedTaxaDict[48] = { "Erinaceus" };
			expextedTaxaDict[44] = { "Oryctolagus" };
			expextedTaxaDict[43] = { "Ochotona" };
			expextedTaxaDict[41] = { "Ictidomys" };
			expextedTaxaDict[37] = { "Rattus" };
			expextedTaxaDict[36] = { "Mus" };
			expextedTaxaDict[35] = { "Dipodomys" };
			expextedTaxaDict[34] = { "Cavia" };
			expextedTaxaDict[32] = { "Tupaia" };
			expextedTaxaDict[29] = { "Tarsius" };
			expextedTaxaDict[27] = { "Otolemur" };
			expextedTaxaDict[26] = { "Microcebus" };
			expextedTaxaDict[23] = { "Macaca" };
			expextedTaxaDict[21] = { "Nomascus" };
			expextedTaxaDict[19] = { "Pongo" };
			expextedTaxaDict[16] = { "Pan" };
			expextedTaxaDict[15] = { "Homo" };
			expextedTaxaDict[14] = { "Gorilla" };
			expextedTaxaDict[13] = { "Callithrix" };
			expextedTaxaDict[11] = { "Dasypus" };
			expextedTaxaDict[10] = { "Choloepus" };
			expextedTaxaDict[7] = { "Procavia" };
			expextedTaxaDict[6] = { "Loxodonta" };
			expextedTaxaDict[5] = { "Echinops" };
			expextedTaxaDict[3] = { "Sarcophilus" };
			expextedTaxaDict[1] = { "Monodelphis" };
			expextedTaxaDict[0] = { "Macropus" };
			vector<int> leafIds = { 77, 65, 64, 63, 62, 60, 59, 58, 56, 54, 52, 50, 49, 48, 44, 43, 41, 37, 36, 35, 34, 32, 29, 27, 26, 23, 21, 19, 16, 15, 
				14, 13, 11, 10, 7 ,6 ,5 ,3 ,1 ,0 };
			for (int i = 0; i < leafIds.size(); ++i)
			{
				expextedChildrenIdsDict[leafIds[i]] = {};
			}

			Assert::AreEqual(expextedChildrenIdsDict.size(), res.second.size());
			for (int i = 0; i < expextedChildrenIdsDict.size(); ++i)
			{
				std::set<int> childrenIds;

				for (Node* child : res.second[i]->children)
				{
					childrenIds.insert(child->id);
				}

				Assert::IsTrue(
					std::equal(
						expextedChildrenIdsDict[i].begin(), expextedChildrenIdsDict[i].end(),
						childrenIds.begin(), childrenIds.end()
					) && expextedChildrenIdsDict[i].size() == childrenIds.size()
				);
			}
			for (int i = 0; i < expextedTaxaDict.size(); ++i)
			{
				Assert::IsTrue(expextedTaxaDict[i] == res.second[i]->keys);
			}
		}

		TEST_METHOD(dpos_for_diff_length_case_e)
		{
			vector<string> profileA = {
				"ATAG",
				"-T-G",
				"CG-G"
			};

			vector<string> profileB = {
				"AT-AG",
				"--T-G",
				"-CG-G"
			};

			double res = compute_distance(profileB, profileA, DistanceType::D_POS);
			Assert::AreEqual(res, 0.333, 0.001);
		}

	};
}
