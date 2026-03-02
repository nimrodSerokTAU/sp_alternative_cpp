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
		
		TEST_METHOD(TestMethod1)
		{
			EvoModel evoModel1(-1, -1, "Blosum50");
			vector<EvoModel> models = { evoModel1 };
			Configuration configuration1(models);
			SPScore sp1(evoModel1, evoModel1.matrix_file_name);
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
			Assert::AreEqual(sum, 228.0); // Example assertion
		}
	};
}
