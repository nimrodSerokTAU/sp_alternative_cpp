
#include <iostream>
#include <fstream>
#include <string>
#include <filesystem>
#include "../sp_alternative_cpp/include/multi_msa_service.h"
#include "../sp_alternative_cpp/include/msa.h"
#include "../sp_alternative_cpp/include/dist_labels_stats.h"
#include "../sp_alternative_cpp/include/distance_calc.h"
using namespace std;
namespace fs = std::filesystem;


void multiple_msa_calc_features_and_labels(string directoryName, string outputDirectory) {

    fs::path dir_path(directoryName);
    fs::path output_dir_path(outputDirectory);
    fs::create_directories(output_dir_path);
    string dir_name = dir_path.filename().string();
    vector<string> files;
    for (const auto& f : fs::directory_iterator(dir_path)) {
        if (f.is_regular_file()) {
            files.push_back(f.path().filename().string());
        }
    }

    auto ordered = get_file_names_ordered(files);

    MSA true_msa(dir_name);
    if (!ordered.true_file_name.empty()) {
        true_msa.read_from_fasta(dir_path / ordered.true_file_name);
    }
    else {
        exit(1);
    }

    const int rows_num = true_msa.sequences.size();
    const int cols_num = true_msa.sequences[0].size();
    vector<vector<int>> true_msa_vectors(cols_num, std::vector<int>(rows_num));
    vector<int> true_msa_map;
    fill_d_seq_vectors(true_msa.sequences, true_msa_vectors, true_msa_map, rows_num, cols_num);

	vector<double> dseq_vectors(ordered.other_file_names.size());
	for (int i = 0; i < ordered.other_file_names.size(); i++) {
		string msa_name = ordered.other_file_names[i];
        cout << msa_name << endl;

        MSA inferred_msa(msa_name);
        inferred_msa.read_from_fasta(dir_path / msa_name);
        inferred_msa.order_sequences(true_msa.seq_names);

        double dSeq = compute_eff_d_seq_from_true(inferred_msa.sequences, true_msa_vectors, true_msa_map);
		dseq_vectors[i] = dSeq;
    }
    ofstream out_file(output_dir_path / (dir_name + "_dseq_from_true.csv"));
	out_file << "code,dseq_from_true\n";
    for (int i = 0; i < ordered.other_file_names.size(); i++) {
        out_file << std::setprecision(4)<< ordered.other_file_names[i] << "," << dseq_vectors[i] << "\n";
	}
}



int main(int argc, char* argv[])
{
    string directoryName = argv[1]; 
	string outputDirectoryName = argv[2];
    multiple_msa_calc_features_and_labels(directoryName, outputDirectoryName);
}


