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
#include "../sp_alternative_cpp/include/neighbor_joining.h"
#include "pch.h"                  
#include <memory>
#include <set>
#include <vector>
#include <algorithm>
#include <cmath>
#include <optional>
#include "../sp_alternative_cpp/include/tree_stats.h"
#include "../sp_alternative_cpp/include/msa.h"
#include "../sp_alternative_cpp/include/w_Sop_stats.h"
#include "../sp_alternative_cpp/include/utils.h"
#include <iostream>
using namespace Microsoft::VisualStudio::CppUnitTestFramework;
using namespace std;

namespace spalternativeUnitTests
{
	string blosum50Path = "D:/code/sp_alternative/sp_alternative/input_config_files/Blosum50.txt";
	string blosum62Path = "D:/code/sp_alternative/sp_alternative/input_config_files/Blosum62.txt";
	string nucPath = "D:/code/sp_alternative/sp_alternative/input_config_files/Nucleotides.txt";
	string newick_of_AATF = "((((Macropus:0.051803,Monodelphis:  0.066021): 0.016682,Sarcophilus:0.068964):0.114355,((Echinops:0.104144,(Loxodonta:0.076474,Procavia:0.076193):0.011550):0.013015,((Choloepus:0.056091,Dasypus:0.040600):0.013681,(((((Callithrix:0.032131,((((Gorilla:0.007042,(Homo:0.002445,Pan:0.002450):0.001237):0.003689,Pongo:0.007508):0.002384,Nomascus:0.015696):0.004674,Macaca:0.013752):0.012205):0.029730,((Microcebus:0.037066,Otolemur:0.050935):0.008091,Tarsius:0.064938):0.007305):0.003377,Tupaia:0.090946):0.000699,(((Cavia:0.116826,(Dipodomys:0.080386,(Mus:0.040313,Rattus:0.033872):0.122329):0.013314):0.000298,Ictidomys:0.062932):0.011064,(Ochotona:0.087746,Oryctolagus:0.057769):0.035947):0.002130):0.006931,(Erinaceus:0.094727,(((((Bos:0.062659,Tursiops:0.024374):0.009782,Sus:0.064336):0.006134,Vicugna:0.049961):0.021643,Equus:0.046559):0.002728,(Sorex:0.126677,((Myotis:0.050452,Pteropus:0.047740):0.006495,(Felis:0.042414,(Canis:0.036675,(Mustela:0.027691,Ailuropoda:0.037553):0.004251):0.008141):0.019485):0.000485):0.003533):0.000528):0.005645):0.012900):0.002853):0.133251):0.019558,Ornithorhynchus:0.195576);";

	UnrootedTree create_unrooted_tree_for_test() {
		vector<unique_ptr<Node>> all_nodes;

		auto n_a = make_unique<Node>(0, set<string>{"a"}, vector<int>{}, 0.0, 0.2);
		auto n_b = make_unique<Node>(1, set<string>{"b"}, vector<int>{}, 0.0, 0.1);
		auto n_c = make_unique<Node>(2, set<string>{"c"}, vector<int>{}, 0.0, 0.15);
		auto n_d = make_unique<Node>(3, set<string>{"d"}, vector<int>{}, 0.0, 0.25);
		auto n_e = make_unique<Node>(4, set<string>{"e"}, vector<int>{}, 0.0, 0.05);

		// Keep raw pointers
		Node* p_a = n_a.get();
		Node* p_b = n_b.get();
		Node* p_c = n_c.get();
		Node* p_d = n_d.get();
		Node* p_e = n_e.get();

		all_nodes.push_back(std::move(n_a));
		all_nodes.push_back(std::move(n_b));
		all_nodes.push_back(std::move(n_c));
		all_nodes.push_back(std::move(n_d));
		all_nodes.push_back(std::move(n_e));

		// --- 2. Internal nodes ---

		// n_a_c
		auto n_a_c = make_unique<Node>(5, set<string>{}, vector<int>{p_a->id, p_c->id}, 0, 0.6);
		Node* p_a_c = n_a_c.get();
		p_a->set_a_father(p_a_c->id);
		p_c->set_a_father(p_a_c->id);
		p_a_c->keys.insert(p_a->keys.begin(), p_a->keys.end());
		p_a_c->keys.insert(p_c->keys.begin(), p_c->keys.end());
		all_nodes.push_back(std::move(n_a_c));


		// n_a_c_d
		auto n_a_c_d = make_unique<Node>(6, set<string>{}, vector<int>{p_a_c->id, p_d->id}, 0, 0.3);
		Node* p_a_c_d = n_a_c_d.get();
		p_a_c->set_a_father(p_a_c_d->id);
		p_d->set_a_father(p_a_c_d->id);
		p_a_c_d->keys.insert(p_a_c->keys.begin(), p_a_c->keys.end());
		p_a_c_d->keys.insert(p_d->keys.begin(), p_d->keys.end());

		all_nodes.push_back(std::move(n_a_c_d));


		// anchor
		auto anchor_node = make_unique<Node>(7, set<string>{}, vector<int>{p_a_c_d->id, p_b->id, p_e->id}, 0);
		Node* anchor = anchor_node.get();

		p_a_c_d->set_a_father(anchor->id);
		p_b->set_a_father(anchor->id);
		p_e->set_a_father(anchor->id);

		for (auto* c : vector<Node*>{ p_a_c_d , p_b, p_e }) {
			anchor->keys.insert(c->keys.begin(), c->keys.end());
		}

		all_nodes.push_back(std::move(anchor_node));
		UnrootedTree tree(anchor, std::move(all_nodes));
		return tree;
	}

	Node* find_node_by_keys(const std::vector<Node*>& nodes, const std::set<std::string>& target_keys) {
		auto it = std::find_if(nodes.begin(), nodes.end(),
			[&](Node* node) {
				return node->keys == target_keys;
			});

		return (it != nodes.end()) ? *it : nullptr;
	}

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
			vector<vector<int>> substitutions_matrix;	
			SPScore::SpSAndGe res = sp1.compute_sp_s_and_sp_ge(profile1, substitutions_matrix, false);
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
			double efficientRes = sp1.compute_efficient_sp(profile1, false);
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
			double efficientRes = sp1.compute_efficient_sp(profile1, false);
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

			double efficientRes = sp1.compute_efficient_sp(profile1, false);
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

			double d_seq = compute_distance(profile1, profile2, DistanceType::D_SEQ);
			Assert::AreEqual(d_seq, 0.611, 0.001);

			double eff_d_seq = compute_eff_d_seq(profile1, profile2);
			Assert::AreEqual(eff_d_seq, 0.611, 0.001);
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
			Assert::AreEqual(res_1_3, 0.295, 0.001);
			double ef_res_t_1 = compute_eff_d_seq(profile1, profile3);
			Assert::AreEqual(res_1_3, 0.295, 0.001);
		}

		TEST_METHOD(dseq_for_diff_methods)
		{
			vector<string> profile1 = {
				"GCATCATT--GT",
				"GC---ATTA-GT",
				"GC---AT-AGGT",
				"G---CAT-AGGT"
			};

			vector<string> profile2 = {
				"GCATCATT-G-T",
				"GCA---TTAG-T",
				"GCA----TAGGT",
				"G---CAT-AGGT"
			};


			double res_1_2 = compute_distance(profile1, profile2, DistanceType::D_SEQ);
			Assert::AreEqual(res_1_2, 0.294, 0.001);
			double ef_res_t_1 = compute_eff_d_seq(profile1, profile2);
			Assert::AreEqual(ef_res_t_1, 0.294, 0.001);
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
			double ef_res_t_1 = compute_eff_d_seq(trueProfile, profileA);
			Assert::AreEqual(ef_res_t_1, 0.091, 0.001);
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
			UnrootedTree ut = UnrootedTree(newick_of_AATF);
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

			vector<Node*> raw_res = get_raw_pointers_from_unique(ut.all_nodes);

			Assert::AreEqual(expextedChildrenIdsDict.size(), raw_res.size());
			for (int i = 0; i < expextedChildrenIdsDict.size(); ++i)
			{
				std::set<int> childrenIds;

				for (int child_id : raw_res[i]->children_ids)
				{
					childrenIds.insert(child_id);
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
				Assert::IsTrue(expextedTaxaDict[i] == raw_res[i]->keys);
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

		TEST_METHOD(neighbor_joining)
		{
			// Distance matrix
			vector<vector<double>> matrix_case_nj = {
				{0.0, 5.0, 9.0, 9.0, 8.0},
				{5.0, 0.0, 10.0, 10.0, 9.0},
				{9.0, 10.0, 0.0, 8.0, 7.0},
				{9.0, 10.0, 8.0, 0.0, 3.0},
				{8.0, 9.0, 7.0, 3.0, 0.0},
			};

			// Keys for leaves
			vector<string> keys_case_nj = { "a", "b", "c", "d", "e" };

			// Create nodes as unique_ptr<Node>
			vector<unique_ptr<Node>> nodesP;
			for (size_t i = 0; i < keys_case_nj.size(); ++i)
			{
				set<string> thisKeySet = { keys_case_nj[i] };
				vector<int> children_ids = {}; // leaf node
				auto new_node = make_unique<Node>(
					static_cast<int>(i),
					thisKeySet,         
					children_ids,       
					1.0,                
					0.0                 
				);

				nodesP.push_back(std::move(new_node)); // MOVE is required!
			}

			// Create NeighborJoining instance (constructor moves nodesP into all_nodes)
			NeighborJoining neighborJoining(matrix_case_nj, std::move(nodesP));

			// Access the calculated tree (tree_res is std::optional)
			auto& tree = neighborJoining.tree_res.value();

			// Get branch lengths and sort
			vector<double> bl_list = tree.get_branches_lengths_list();
			sort(bl_list.begin(), bl_list.end());

			// Expected branch lengths
			vector<double> expected_bl_list = { 1.0, 2.0, 2.0, 2.0, 3.0, 3.0, 4.0 };

			// Check sizes
			Assert::AreEqual(expected_bl_list.size(), bl_list.size());

			// Compare each branch length (allow small epsilon for floating-point)
			const double EPS = 1e-6;
			for (size_t i = 0; i < bl_list.size(); ++i)
			{
				Assert::IsTrue(fabs(expected_bl_list[i] - bl_list[i]) < EPS);
			}
			string newick_str = tree.print_newick();
			Assert::AreEqual(newick_str.c_str(), "(((a:2.000000,b:3.000000):3.000000,c:4.000000):2.000000,d:2.000000,e:1.000000):0.000000;");
		}

		TEST_METHOD(neighbor_joining_case_b)
		{
			// Distance matrix
			vector<vector<double>> matrix_case_nj = {
				{0.000, 0.385, 0.385, 0.385, 0.692, 0.615, 0.769, 0.538, 0.615},
				{0.385, 0.000, 0.231, 0.000, 0.538, 0.462, 0.692, 0.385, 0.538},
				{0.385, 0.231, 0.000, 0.231, 0.308, 0.231, 0.538, 0.231, 0.308},
				{0.385, 0.000, 0.231, 0.000, 0.538, 0.462, 0.692, 0.385, 0.538},
				{0.692, 0.538, 0.308, 0.538, 0.000, 0.385, 0.231, 0.462, 0.462},
				{0.615, 0.462, 0.231, 0.462, 0.385, 0.000, 0.615, 0.385, 0.077},
				{0.769, 0.692, 0.538, 0.692, 0.231, 0.615, 0.000, 0.462, 0.615},
				{0.538, 0.385, 0.231, 0.385, 0.462, 0.385, 0.462, 0.000, 0.462},
				{0.615, 0.538, 0.308, 0.538, 0.462, 0.077, 0.615, 0.462, 0.000}
			};

			vector<string> keys_case_nj = { "rayfinfish", "frog", "turtle", "salamander", "crocodile", "lizard", "bird", "mammal", "snake"};

			vector<unique_ptr<Node>> nodesP;
			for (size_t i = 0; i < keys_case_nj.size(); ++i)
			{
				set<string> thisKeySet = { keys_case_nj[i] };
				vector<int> children_ids = {}; // leaf node
				auto new_node = make_unique<Node>(
					static_cast<int>(i),
					thisKeySet,
					children_ids,
					1.0,
					0.0
				);

				nodesP.push_back(std::move(new_node));
			}

			NeighborJoining neighborJoining(matrix_case_nj, std::move(nodesP));
			auto& tree = neighborJoining.tree_res.value();
			string newick_str = tree.print_newick();
			Assert::AreEqual(newick_str.c_str(), "(((((frog:0.000000,salamander:0.000000):0.125313,rayfinfish:0.259687):0.089542,turtle:0.025958):0.038094,mammal:0.154156):0.028969,(bird:0.186786,crocodile:0.044214):0.187344,(lizard:0.009792,snake:0.067208):0.177906):0.000000;");
		}

		TEST_METHOD(parsimony)
		{
			vector<unique_ptr<Node>> all_nodes;

			auto n_a = make_unique<Node>(0, set<string>{"a"}, vector<int>{}, 0.0);
			auto n_b = make_unique<Node>(1, set<string>{"b"}, vector<int>{}, 0.0);
			auto n_c = make_unique<Node>(2, set<string>{"c"}, vector<int>{}, 0.0);
			auto n_d = make_unique<Node>(3, set<string>{"d"}, vector<int>{}, 0.0);
			auto n_e = make_unique<Node>(4, set<string>{"e"}, vector<int>{}, 0.0);

			// Keep raw pointers
			Node* p_a = n_a.get();
			Node* p_b = n_b.get();
			Node* p_c = n_c.get();
			Node* p_d = n_d.get();
			Node* p_e = n_e.get();

			all_nodes.push_back(std::move(n_a));
			all_nodes.push_back(std::move(n_b));
			all_nodes.push_back(std::move(n_c));
			all_nodes.push_back(std::move(n_d));
			all_nodes.push_back(std::move(n_e));

			// --- 2. Internal nodes ---

			// n_a_b
			auto n_a_b = make_unique<Node>(5, set<string>{}, vector<int>{p_a->id, p_b->id}, 0);
			Node* p_a_b = n_a_b.get();
			p_a->set_a_father(p_a_b->id);
			p_b->set_a_father(p_a_b->id);
			p_a_b->keys.insert(p_a->keys.begin(), p_a->keys.end());
			p_a_b->keys.insert(p_b->keys.begin(), p_b->keys.end());

			all_nodes.push_back(std::move(n_a_b));

			// n_a_b_c
			auto n_a_b_c = make_unique<Node>(6, set<string>{}, vector<int>{p_a_b->id, p_c->id}, 0);
			Node* p_a_b_c = n_a_b_c.get();
			p_a_b->set_a_father(p_a_b_c->id);
			p_c->set_a_father(p_a_b_c->id);
			p_a_b_c->keys.insert(p_a_b->keys.begin(), p_a_b->keys.end());
			p_a_b_c->keys.insert(p_c->keys.begin(), p_c->keys.end());

			all_nodes.push_back(std::move(n_a_b_c));

			// anchor
			auto anchor_node = make_unique<Node>(7, set<string>{}, vector<int>{p_a_b_c->id, p_d->id, p_e->id}, 0);
			Node* anchor = anchor_node.get();

			p_a_b_c->set_a_father(anchor->id);
			p_d->set_a_father(anchor->id);
			p_e->set_a_father(anchor->id);

			for (auto* c : vector<Node*>{ p_a_b_c , p_d, p_e }) {
				anchor->keys.insert(c->keys.begin(), c->keys.end());
			}

			all_nodes.push_back(std::move(anchor_node));

			// --- 3. Alignment ---
			vector<string> aln = {
				"AYCDDDW",
				"AVVDDDW",
				"AYCDDDW",
				"AVVDDDW",
				"APVDDDW"
			};

			vector<string> names = { "a", "b", "c", "d", "e" };

			// --- 4. Build tree ---
			UnrootedTree tree(anchor, std::move(all_nodes));

			// --- 5. Run parsimony ---
			vector<int> res = calc_parsimony(tree, aln, names);

			// --- 6. Expected result ---
			vector<int> expected = { 0, 3, 2, 0, 0, 0, 0 };

			// --- 7. Assert ---
			Assert::AreEqual(res.size(), expected.size());
			for (size_t i = 0; i < res.size(); i++) {
				Assert::AreEqual(res[i], expected[i]);
			}
		}

		TEST_METHOD(henikoff_w)
		{
			auto msa_ptr = make_unique<MSA>("test");
			msa_ptr->add_sequence(string{"AT-CGC"}, "a");
			msa_ptr->add_sequence(string{"ACATG-"}, "b");
			msa_ptr->add_sequence(string{"AT-CG-"}, "c");
			msa_ptr->add_sequence(string{"ATC-GA"}, "d");
			msa_ptr->add_sequence(string{"TTATGC"}, "e");
			auto w_sop_ptr = make_unique<WSopStats>(msa_ptr->dataset_name, msa_ptr->get_taxa_num(), msa_ptr->get_msa_len());
			std::pair<std::vector<double>, std::vector<double>> seq_weights = w_sop_ptr->compute_seq_w_henikoff_vars(msa_ptr->sequences);
			vector<double> seq_weights_no_gap_expected = {
				0.15454545454545454,
				0.22272727272727275,
				0.10909090909090909,
				0.24545454545454548,
				0.2681818181818182,
			};
			vector<double> seq_weights_with_gap_expected = {
				0.15714285714285717,
				0.21071428571428572,
				0.15714285714285717,
				0.2642857142857143,
				0.21071428571428572,
			};
				
			for (size_t i = 0; i < seq_weights_with_gap_expected.size(); ++i)
			{
				Assert::AreEqual(seq_weights.first[i], seq_weights_with_gap_expected[i]);
			}

			for (size_t i = 0; i < seq_weights_no_gap_expected.size(); ++i)
			{
				Assert::AreEqual(seq_weights.second[i], seq_weights_no_gap_expected[i]);
			}
		}

		TEST_METHOD(henikoff_with_gaps_value)
		{
			auto msa_ptr = make_unique<MSA>("test");
			msa_ptr->add_sequence(string{ "AT-CGC" }, "a");
			msa_ptr->add_sequence(string{ "ACATG-" }, "b");
			msa_ptr->add_sequence(string{ "AT-CG-" }, "c");
			msa_ptr->add_sequence(string{ "ATC-GA" }, "d");
			msa_ptr->add_sequence(string{ "TTATGC" }, "e");
			
			EvoModel evoModel1(-10, -0.5, "Blosum50");
			vector<EvoModel> models = { evoModel1 };
			Configuration config(models, SopCalcTypes::EFFICIENT, "", "", "", { WeightMethods::HENIKOFF_WG });
			auto w_sop_ptr = make_unique<WSopStats>(msa_ptr->dataset_name, msa_ptr->get_taxa_num(), msa_ptr->get_msa_len());
			w_sop_ptr->calc_seq_weights(config.additional_weights, msa_ptr->sequences, msa_ptr->seq_names, UnrootedTree());
			SPScore sp(evoModel1, blosum62Path);
			w_sop_ptr->calc_w_sp(msa_ptr->sequences, sp);
			double henikoff_with_gaps = w_sop_ptr->sp_HENIKOFF_with_gaps;
			double expected_value = -1.682110969387752;
			Assert::IsTrue(abs(henikoff_with_gaps - expected_value) < 1e-10);
		}

		TEST_METHOD(clustalw_for_a_tree)
		{
			auto msa_ptr = make_unique<MSA>("test");
			msa_ptr->add_sequence(string{ "AT-CGC" }, "a");
			msa_ptr->add_sequence(string{ "ACATG-" }, "b");
			msa_ptr->add_sequence(string{ "AT-CG-" }, "c");
			msa_ptr->add_sequence(string{ "ATC-GA" }, "d");
			msa_ptr->add_sequence(string{ "TTATGC" }, "e");

			string newick = "(((a:1,b:2):3,c:4):5,d:6,e:7);";
			EvoModel evoModel1(-10, -0.5, "Blosum62");
			vector<EvoModel> models = { evoModel1 };
			Configuration config(models, SopCalcTypes::EFFICIENT, "", "", "", { WeightMethods::CLUSTAL_MID_ROOT });

			UnrootedTree ut = UnrootedTree(newick);
			auto w_sop_ptr = make_unique<WSopStats>(msa_ptr->dataset_name, msa_ptr->get_taxa_num(), msa_ptr->get_msa_len());
			w_sop_ptr->calc_seq_weights(config.additional_weights, msa_ptr->sequences, msa_ptr->seq_names, ut);
			SPScore sp(evoModel1, blosum62Path);
			w_sop_ptr->calc_w_sp(msa_ptr->sequences, sp);
			double sp_clustul_mid_root = w_sop_ptr->sp_CLUSTAL_WEIGHTS_mid_root;
			double expected_value = -1.5516;
			Assert::IsTrue(abs(sp_clustul_mid_root - expected_value) < 1e-4);
		}

		TEST_METHOD(calc_dseq_from_file)
		{
			string trueProfilePath = "true_msa.fas";
			string profileAPath = "test_msa.fas";
			
			MSA true_msa("trueMSA");
			true_msa.read_from_fasta(trueProfilePath);

			MSA test_msa("testMSA");
			test_msa.read_from_fasta(profileAPath);

			double ef_res_t_1 = compute_eff_d_seq(true_msa.sequences, test_msa.sequences);
			Assert::AreEqual(ef_res_t_1, 0.044, 0.001);
			
			double res_t_1 = compute_distance(true_msa.sequences, test_msa.sequences, DistanceType::D_SEQ);
			Assert::AreEqual(res_t_1, 0.044, 0.001);
		}

		TEST_METHOD(rooting_case_a)
		{
			UnrootedTree unrooted = create_unrooted_tree_for_test();
			RootingPoint rp(6, 7, 0.22, 0.08, -1.0);
			RootedTree tree = RootedTree(unrooted, rp);

			vector<double> expected_branch_lengths = { 0.2, 0.1, 0.15, 0.25, 0.05, 0.6, 0.22, 0.08, 0 };
			for (int i = 0; i < 9; ++i) {
				Node* thisNode = tree.all_nodes[i].get();
				Assert::AreEqual(thisNode->branch_length, expected_branch_lengths[i], 0.001);
				if (i < 5) {
					Assert::AreEqual(int(thisNode->children_ids.size()), 0);
				}
				else {
					Assert::AreEqual(int(thisNode->children_ids.size()), 2);
				}
			}
		}

		TEST_METHOD(rooting_case_b)
		{
			UnrootedTree unrooted = create_unrooted_tree_for_test();
			RootingPoint rp(6, 5, 0.2, 0.4, -1.0);
			RootedTree tree = RootedTree(unrooted, rp);

			vector<double> expected_branch_lengths = { 0.2, 0.1, 0.15, 0.25, 0.05, 0.4, 0.2, 0.3, 0 };
			for (int i = 0; i < 9; ++i) {
				Node* thisNode = tree.all_nodes[i].get();
				Assert::AreEqual(thisNode->branch_length, expected_branch_lengths[i], 0.001);
				if (i < 5) {
					Assert::AreEqual(int(thisNode->children_ids.size()), 0);
				}
				else {
					Assert::AreEqual(int(thisNode->children_ids.size()), 2);
				}
			}
		}
		
		TEST_METHOD(rooting_case_c)
		{
			UnrootedTree unrooted = create_unrooted_tree_for_test();
			RootingPoint rp(0, 5, 0.07, 0.13, -1.0);
			RootedTree tree = RootedTree(unrooted, rp);

			vector<double> expected_branch_lengths = { 0.07, 0.1, 0.15, 0.25, 0.05, 0.13, 0.6, 0.3, 0 };
			for (int i = 0; i < 9; ++i) {
				Node* thisNode = tree.all_nodes[i].get();
				Assert::AreEqual(thisNode->branch_length, expected_branch_lengths[i], 0.001);
				if (i < 5) {
					Assert::AreEqual(int(thisNode->children_ids.size()), 0);
				}
				else {
					Assert::AreEqual(int(thisNode->children_ids.size()), 2);
				}
			}
		}

		TEST_METHOD(rooting_case_d)
		{
			//string newick = "(((((((((((((((Grammomys_surdaster:0.094018,Rattus_norvegicus:0.091931):0.091668,Apodemus_sylvaticus:0.186518):0.130051,(Mesocricetus_auratus:0.174981,Phodopus_roborovskii:0.178931):0.190541):0.054003,Cricetulus_griseus:0.409168):0.127290,(((Mandrillus_leucophaeus:0.403218,Mus_pahari:0.414545):0.074584,Tupaia_chinensis:0.416321):0.044631,Hipposideros_armiger:0.503176):0.007727):0.027741,((((Dipodomys_spectabilis:0.597143,Perognathus_longimembris_pacificus:1.493979):0.414428,Jaculus_jaculus:0.540011):0.495539,Equus_przewalskii:0.504461):0.042163,Manis_pentadactyla:0.488968):0.050030):0.009218,Myotis_davidii:0.532629):0.059763,Condylura_cristata:0.430158):0.025425,Equus_quagga:0.472891):0.063130,(Balaenoptera_acutorostrata_scammoni:0.275313,Balaenoptera_musculus:0.213649):0.033392):0.043378,Hyaena_hyaena:0.309704):0.054654,(((((((((Rhinolophus_ferrumequinum:0.099350,Sus_scrofa:0.133844):0.020719,Prionailurus_bengalensis:0.168304):0.014171,(Carlito_syrichta:0.154886,Otolemur_garnettii:0.144293):0.052677):0.018333,Suricata_suricatta:0.212266):0.022806,Dasypus_novemcinctus:0.242557):0.068824,(Galeopterus_variegatus:0.122153,Macaca_nemestrina:0.156135):0.166575):0.009390,(Aotus_nancymaae:0.297897,Cebus_imitator:0.336023):0.075308):0.023396,Microcebus_murinus:0.373571):0.044918,Rousettus_aegyptiacus:0.403349):0.029021):0.058491,Molossus_molossus:0.435736):0.028099,(((Myodes_glareolus:0.072444,Peromyscus_leucopus:0.060376):0.179877,Sorex_araneus:0.337851):0.106576,(Acomys_russatus:0.281378,Nannospalax_galili:0.248545):0.109855):0.068974):0.015453,((((((((((((((Sturnira_hondurensis:0.944934,Urocitellus_parryii:1.055066):0.055960,Erinaceus_europaeus:0.944040):0.271553,Elephantulus_edwardii:0.728447):0.175118,Vulpes_vulpes:0.881384):0.164137,Vicugna_pacos:0.663009):0.088956,(Ceratotherium_simum_simum:0.644495,Loxodonta_africana:0.716172):0.178792):0.059907,(Cercocebus_atys:0.919082,Suncus_etruscus:0.726776):0.032272):0.026291,(((Choloepus_didactylus:0.476678,Chrysochloris_asiatica:0.665470):0.239322,Trichechus_manatus_latirostris:0.611229):0.052494,(Ochotona_princeps:0.711918,Odocoileus_virginianus_texanus:0.780655):0.085818):0.053783):0.055438,(((Desmodus_rotundus:0.287195,Phyllostomus_discolor:0.370743):0.081775,Artibeus_jamaicensis:0.602598):0.057247,Bos_indicus:0.512644):0.027219):0.024561,Miniopterus_natalensis:0.701458):0.053564,Castor_canadensis:0.665204):0.007626,(Orycteropus_afer_afer:0.559711,Phascolarctos_cinereus:0.639816):0.035858):0.057740,Talpa_occidentalis:0.534694):0.046753,((((Meles_meles:0.337363,Mirounga_leonina:0.158289):0.120636,Vulpes_lagopus:0.282327):0.129021,Nomascus_leucogenys:0.359029):0.189430,Ornithorhynchus_anatinus:0.560461):0.052013):0.013405,Physeter_catodon:0.440763);";
			string newick = "(((((((((Dasypus_novemcinctus:0.242557,(((Carlito_syrichta:0.154886,Otolemur_garnettii:0.144293):0.052677,((Rhinolophus_ferrumequinum:0.099350,Sus_scrofa:0.133844):0.020719,Prionailurus_bengalensis:0.168304):0.014171):0.018333,Suricata_suricatta:0.212266):0.022806):0.068824,(Macaca_nemestrina:0.156135,Galeopterus_variegatus:0.122153):0.166575):0.009390,(Aotus_nancymaae:0.297897,Cebus_imitator:0.336023):0.075308):0.023396,Microcebus_murinus:0.373571):0.044918,Rousettus_aegyptiacus:0.403349):0.029021,((((((((((Perognathus_longimembris_pacificus:1.493979,Dipodomys_spectabilis:0.597143):0.414428,Jaculus_jaculus:0.540011):0.495539,Equus_przewalskii:0.504461):0.042163,Manis_pentadactyla:0.488968):0.050030,((((Apodemus_sylvaticus:0.186518,(Grammomys_surdaster:0.094018,Rattus_norvegicus:0.091931):0.091668):0.130051,(Phodopus_roborovskii:0.178931,Mesocricetus_auratus:0.174981):0.190541):0.054003,Cricetulus_griseus:0.409168):0.127290,(((Mus_pahari:0.414545,Mandrillus_leucophaeus:0.403218):0.074584,Tupaia_chinensis:0.416321):0.044631,Hipposideros_armiger:0.503176):0.007727):0.027741):0.009218,Myotis_davidii:0.532629):0.059763,Condylura_cristata:0.430158):0.025425,Equus_quagga:0.472891):0.063130,(Balaenoptera_musculus:0.213649,Balaenoptera_acutorostrata_scammoni:0.275313):0.033392):0.043378,Hyaena_hyaena:0.309704):0.054654):0.058491,Molossus_molossus:0.435736):0.028099,((Acomys_russatus:0.281378,Nannospalax_galili:0.248545):0.109855,((Peromyscus_leucopus:0.060376,Myodes_glareolus:0.072444):0.179877,Sorex_araneus:0.337851):0.106576):0.068974):0.015453,((((Castor_canadensis:0.665204,(((((Ochotona_princeps:0.711918,Odocoileus_virginianus_texanus:0.780655):0.085818,((Chrysochloris_asiatica:0.665470,Choloepus_didactylus:0.476678):0.239322,Trichechus_manatus_latirostris:0.611229):0.052494):0.053783,(((((((Urocitellus_parryii:1.055066,Sturnira_hondurensis:0.944934):0.055960,Erinaceus_europaeus:0.944040):0.271553,Elephantulus_edwardii:0.728447):0.175118,Vulpes_vulpes:0.881384):0.164137,Vicugna_pacos:0.663009):0.088956,(Ceratotherium_simum_simum:0.644495,Loxodonta_africana:0.716172):0.178792):0.059907,(Cercocebus_atys:0.919082,Suncus_etruscus:0.726776):0.032272):0.026291):0.055438,((Artibeus_jamaicensis:0.602598,(Phyllostomus_discolor:0.370743,Desmodus_rotundus:0.287195):0.081775):0.057247,Bos_indicus:0.512644):0.027219):0.024561,Miniopterus_natalensis:0.701458):0.053564):0.007626,(Orycteropus_afer_afer:0.559711,Phascolarctos_cinereus:0.639816):0.035858):0.057740,Talpa_occidentalis:0.534694):0.046753,((Nomascus_leucogenys:0.359029,((Mirounga_leonina:0.158289,Meles_meles:0.337363):0.120636,Vulpes_lagopus:0.282327):0.129021):0.189430,Ornithorhynchus_anatinus:0.560461):0.052013):0.013405,Physeter_catodon:0.440763);";
			UnrootedTree ut = UnrootedTree(newick);
			vector<Node*> raw_nodes = get_raw_pointers_from_unique(ut.all_nodes);
			Node* startNode = find_node_by_keys(raw_nodes, { "Tupaia_chinensis", "Rattus_norvegicus", "Dipodomys_spectabilis", "Perognathus_longimembris_pacificus", "Mesocricetus_auratus", "Grammomys_surdaster", "Mus_pahari", "Hipposideros_armiger", "Equus_przewalskii", "Cricetulus_griseus", "Phodopus_roborovskii", "Mandrillus_leucophaeus", "Apodemus_sylvaticus", "Jaculus_jaculus", "Manis_pentadactyla" });
			Node* endNode = find_node_by_keys(raw_nodes, { "Myotis_davidii", "Tupaia_chinensis", "Rattus_norvegicus", "Dipodomys_spectabilis", "Perognathus_longimembris_pacificus", "Mesocricetus_auratus", "Grammomys_surdaster", "Mus_pahari", "Hipposideros_armiger", "Equus_przewalskii", "Cricetulus_griseus", "Phodopus_roborovskii", "Mandrillus_leucophaeus", "Apodemus_sylvaticus", "Jaculus_jaculus", "Manis_pentadactyla" });

			RootingPoint rp(startNode->id, endNode->id, 0.008774422290135586, 0.00044397946182012094, -1.0);
			RootedTree tree = RootedTree(ut, rp);

			vector<double> expected_branch_lengths = { 0.242556,1.493979,0.597142,0.540011,0.281377,0.186517,0.414544,0.094017,0.091930,0.178931,0.174981,0.409167,0.060375,0.072444,0.248545,
				0.665204,0.711918,1.055066,0.359029,0.156134,0.919082,0.403217,0.297897,0.336023,0.154886,0.144293,0.373571,0.122152,0.416320,0.435735,0.532629,0.701457,0.944933,0.602598,
				0.370742,0.287194,0.503176,0.099350,0.403349,0.158289,0.337362,0.282327,0.881383,0.212265,0.168303,0.309703,0.488967,0.440763,0.213649,0.275312,0.512644,0.780655,0.133843,
				0.663008,0.644494,0.504461,0.472890,0.534693,0.430157,0.944039,0.337850,0.726776,0.728447,0.665470,0.611228,0.716171,0.559711,0.639815,0.560461,0.476677,0.414428,0.055960,
				0.271552,0.495538,0.239321,0.175118,0.091668,0.190541,0.130050,0.120636,0.178792,0.164136,0.179876,0.054002,0.129021,0.074584,0.127290,0.042162,0.044631,0.050029,0.088955,
				0.007726,0.027741,0.189430,0.008774,0.000443,0.106576,0.059763,0.109854,0.025425,0.166574,0.085817,0.059907,0.052493,0.032271,0.026291,0.081775,0.053783,0.052677,0.020718,
				0.014171,0.018333,0.022805,0.075307,0.068824,0.033392,0.057246,0.055438,0.063130,0.068973,0.052013,0.009390,0.023395,0.043378,0.044917,0.027219,0.024560,0.029021,0.053563,
				0.054653,0.035858,0.007625,0.057740,0.046753,0.058490,0.028099,0.013404,0.015453,0 };
			for (int i = 0; i < 139; ++i) {
				if (i == 90) {
					bool debug = true;
				}
				Node* thisNode = tree.all_nodes[i].get();
				cout << "Node " << i << ": branch_length = " << thisNode->branch_length << ", expected = " << expected_branch_lengths[i] << endl;
				//Assert::AreEqual(thisNode->branch_length, expected_branch_lengths[i], 0.001);
				int children_ids_size = int(thisNode->children_ids.size());
				Assert::IsTrue(children_ids_size == 0 || children_ids_size == 2);
			}

			set<string> root_keys = tree.get_my_keys_set(tree.root);
			Assert::IsTrue(root_keys.size() == 70);
			int root_children_count = tree.get_my_children_count(tree.root);
			Assert::IsTrue(root_children_count == 138);

			vector<Node*> tree_raw_nodes = get_raw_pointers_from_unique(tree.all_nodes);


			int child_a_id = tree.root->children_ids[0];
			int child_b_id = tree.root->children_ids[1];

			set<string> c_a_keys = tree.get_my_keys_set(tree_raw_nodes[child_a_id]);
			Assert::IsTrue(c_a_keys.size() == 15);
			set<string> c_b_keys = tree.get_my_keys_set(tree_raw_nodes[child_b_id]);
			Assert::IsTrue(c_b_keys.size() == 55);

			string newick_str = tree.print_newick();
			bool a = 1;
		}

		TEST_METHOD(rooting_case_e)
		{
			string msaFilePath = "bali_phy_msa.1.fas";
			string dataset_name = "test";
			MSA msa(dataset_name);
			msa.read_from_fasta(msaFilePath);
			msa.build_nj_tree();
			UnrootedTree* ut = msa.tree.get();

			vector<Node*> raw_nodes;
			raw_nodes.reserve(ut->all_nodes.size());

			transform(ut->all_nodes.begin(), ut->all_nodes.end(),
				back_inserter(raw_nodes),
				[](const unique_ptr<Node>& p) { return p.get(); });


			RootingPoint rp = get_rooting_point(RootingMethods::LONGEST_PATH_MID, *ut);
			RootedTree rt = RootedTree(*ut, rp);
			string x = rt.print_newick();

			vector<Node*> rt_raw_nodes = get_raw_pointers_from_unique(rt.all_nodes);
			
			for (int i = 0; i < 139; ++i) {
				Node* thisNode = rt_raw_nodes[i];
				//Assert::AreEqual(thisNode->branch_length, expected_branch_lengths[i], 0.001);
				int children_ids_size = int(thisNode->children_ids.size());
				Assert::IsTrue(children_ids_size == 0 || children_ids_size == 2);
			}

			set<string> root_keys = rt.get_my_keys_set(rt.root);
			Assert::IsTrue(root_keys.size() == 70);
			int root_children_count = rt.get_my_children_count(rt.root);
			Assert::IsTrue(root_children_count == 138);

			int child_a_id = rt.root->children_ids[0];
			int child_b_id = rt.root->children_ids[1];

			set<string> c_a_keys = rt.get_my_keys_set(rt_raw_nodes[child_a_id]);
			Assert::IsTrue(c_a_keys.size() == 15);
			set<string> c_b_keys = rt.get_my_keys_set(rt_raw_nodes[child_b_id]);
			Assert::IsTrue(c_b_keys.size() == 55);

			string newick_str = rt.print_newick();
			bool a = 1;
		}

		TEST_METHOD(compute_sp_s_and_sp_ge_nucleotides)
		{
			EvoModel evoModel1(-3, -1, "Nucleotides");
			vector<EvoModel> models = { evoModel1 };
			Configuration configuration1(models);
			configuration1.is_unified_file = true;
			SPScore sp1(evoModel1, nucPath);
			vector<string> profile1 = {
				"ATGTC---GT",
				"AA-TCG--AT",
				"AA--CGAGGT"
			};
			vector<vector<int>> substitutions_matrix;
			SPScore::EfficientSpParts res = sp1.compute_efficient_sp_parts(profile1, substitutions_matrix, true);

			vector<double> vw1 = { 1.0, 1.0, 1.0 };
			vector<vector<double>> vvw1 = { vw1 };
			vector<double> resNaive = sp1.compute_naive_sp_score(profile1, &vvw1);
			double totalScore = res.sp_match_score + res.sp_mismatch_score + res.go_score + res.sp_score_gap_e;
			Assert::AreEqual(totalScore, resNaive[0]);
			Assert::AreEqual(res.sp_mismatch_score + res.sp_match_score, 9.0);

			vector<int> A_vec = { 0, 2, 2, 0 };
			vector<int> T_vec = { 2, 0, 0, 0 };
			vector<int> G_vec = { 2, 0, 0, 0 };
			vector<int> C_vec = { 0, 0, 0, 0 };
			vector<vector<int>> expected_subs_matrix = { A_vec , T_vec, G_vec, C_vec };

			vector<string> nucs = { "A", "T", "G", "C" };

			for (int i = 0; i < nucs.size(); ++i) {
				vector<int> expected_row = expected_subs_matrix[i];
				for (int j = i + 1; j < nucs.size(); ++j) {
					Assert::AreEqual(substitutions_matrix[i][j], expected_row[j]);
				}
			}
		}

		TEST_METHOD(print_nucleotides_res)
		{
			EvoModel evoModel1(-3, -1, "Nucleotides");
			vector<EvoModel> models = { evoModel1 };
			string trueProfilePath = "nuc_TRUE.fas";
			string profileAPath = "nuc_test.fasta";
			
			Configuration configuration(
				models,
				SopCalcTypes::EFFICIENT,
				"",
				"D:/code/sp_alternative/sp_alternative/nuc_output",
				"D:/code/sp_alternative/sp_alternative/input_config_files",
				{},
				{ 8, 16, 32 },
				{ StatsOutput::ALL },
				true
			);

			SPScore sp1(evoModel1, nucPath);
			vector<SPScore> sp_models = { sp1 };
			for (const auto& m : models) {
				fs::path matrix_path = fs::path(configuration.matrix_dir_path) / m.matrix_file_name;
				matrix_path.replace_extension(".txt");
				sp_models.emplace_back(m, matrix_path);
			}
			cout << "End subs matrices: " << endl;

			fs::path output_dir_path = fs::path(configuration.output_file_dir_path);
			fs::create_directories(output_dir_path);

			MSA true_msa("true");
			true_msa.read_from_fasta(trueProfilePath);
			
			MSA inferred_msa("nuc_test");
			inferred_msa.read_from_fasta(profileAPath);
			inferred_msa.order_sequences(true_msa.seq_names);
			inferred_msa.calc_and_print_stats(true_msa, configuration, sp_models, output_dir_path,	true_msa.tree.get(), true);

		}
	};
}
