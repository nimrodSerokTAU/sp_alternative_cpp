#include "pch.h"
#include "CppUnitTest.h"
#include "../sp_alternative_cpp/include/evo_model.h"
#include <vector>
#include "../sp_alternative_cpp/include/config.h"
#include "../sp_alternative_cpp/include/sp_score.h"
#include <numeric>
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
	};
}
