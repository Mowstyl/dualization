#include "bipartite_02.hpp"
#include "bipartite_hypergraph.hpp"
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <mpi.h>
#include <random>
#include <sstream>
#include <string>
#include <vector>

using namespace std;

// Write a row to a CSV file
void write_csv_row(ofstream &file, const vector<string> &row) {
    for (size_t i = 0; i < row.size(); ++i) {
        file << row[i];
        if (i + 1 < row.size())
            file << ";";
    }
    file << "\n";
}

// Implementation of the algorithm execution
vector<string> execute_algorithm(int algorithm_num, bipartite_graph &b,
                                 hypergraph &g, int world_rank, int world_size) {
    string name;
    double elapsed = 0.0;
    string output = "";
    // Prepare data for algorithm_A: convert psi (set of vector) to vector of
    // vector
    vector<vector<int>> F(g.psi.begin(), g.psi.end());

    auto start = chrono::high_resolution_clock::now();
    if (algorithm_num == 0) {
        name = "algorithm A";
        output = to_string(algorithm_A(F, F));
    } else if (algorithm_num == 1) {
        name = "hitting set algorithm";
        output = b.compute_number_hitting_sets().get_str();
    } else if (algorithm_num == 2) {
        name = "search x";
        auto search_result = g.search_x();
        ostringstream oss;
        oss << (search_result.first ? "self-dual" : "not self-dual");
        if (!search_result.first) {
            oss << " [";
            for (size_t i = 0; i < search_result.second.size(); ++i) {
                oss << to_string(search_result.second[i]);
                if (i + 1 < search_result.second.size())
                    oss << ",";
            }
            oss << "]";
        }
        output = oss.str();
    } else if (algorithm_num == 3) {
        name = "search x parallel";
        auto search_result = g.search_x_par(world_rank, world_size);
        ostringstream oss;
        oss << (search_result.first ? "self-dual" : "not self-dual");
        if (!search_result.first) {
            oss << " [";
            for (size_t i = 0; i < search_result.second.size(); ++i) {
                oss << to_string(search_result.second[i]);
                if (i + 1 < search_result.second.size())
                    oss << ",";
            }
            oss << "]";
        }
        output = oss.str();
    } else if (algorithm_num == 4) {
        name = "find_uncovering_hitting_set";
        vector<vector<int>> H = g.get_ordered_list_psi();
        bipartite_02 b2(H);
        bool found = b2.find_uncovering_hitting_set(0);
        ostringstream oss;
        oss << (!found ? "self-dual" : "not self-dual");
        if (found) {
            int i = 0;
            oss << " [";
            for (int v : b2.counter_example) {
                if (i > 0)
                    oss << ",";
                i++;
                oss << to_string(v);
            }
            oss << "]";
        }
        output = oss.str();
    } else if (algorithm_num == 5) {
        name = "sum_f";
        output = g.sum_f().get_str();
    } else {
        if (world_rank == 0)
            cerr << "Unknown algorithm number: " << algorithm_num << endl;
        return {};
    }

    auto end = chrono::high_resolution_clock::now();
    elapsed = chrono::duration<double>(end - start).count();

    if (world_rank == 0)
        cout << name << ": " << elapsed << " " << output << endl;
    return {name, to_string(elapsed), output};
}

void run_multiple_test(const string &fname, const vector<string> &hypergraphs, int world_rank, int world_size) {
    ofstream csvFile;
    if (world_rank == 0) {
        csvFile.open(fname, ios::app);
        if (!csvFile.is_open()) {
            if (world_rank == 0)
                cerr << "Impossibile to open file " << fname << endl;
            return;
        }
    }

    hypergraph g(0);
    for (const string &filepath : hypergraphs) {
        string filename = filesystem::path(filepath).filename().string();
        string datetime_string = get_datetime_string();
        g.load(filepath);
        int n = g.n;
        auto st = g.stat();
        bipartite_graph b(g.get_ordered_list_psi());
        if (world_rank == 0)
            cout << "hypergraph : " << filename << endl;
        if (world_rank == 0)
            cout << "Number of vertices: " << n
                 << "; Number of hyperedges : " << b.Edges.size() << endl;

        vector<string> row_head = {filename,
                                   datetime_string,
                                   to_string(static_cast<int>(st["n"])),
                                   to_string(static_cast<int>(st["m"])),
                                   to_string(static_cast<int>(st["min_edge"])),
                                   to_string(static_cast<int>(st["max_edge"])),
                                   to_string(st["avg"])};

        for (int algorithm_num = 0; algorithm_num < 6; ++algorithm_num) {
            if (algorithm_num != 3 && world_rank != 0)
                continue;
            auto row = execute_algorithm(algorithm_num, b, g, world_rank, world_size);
            vector<string> new_row = row_head;
            new_row.insert(new_row.end(), row.begin(), row.end());
            write_csv_row(csvFile, new_row);
        }

        if (world_rank == 0)
            cout << "----------------------------------------" << endl;
    }
    csvFile.close();
}

bool ends_with(string const &fullString, string const &ending) {
    if (fullString.length() >= ending.length())
        return (0 == fullString.compare(fullString.length() - ending.length(),
                                        ending.length(), ending));
    else
        return false;
}

int main(int argc, char *argv[]) {
    if (argc < 2) {
        cout << "Syntax: " << argv[0]
             << " <hypergraph json 1> ... <hypergraph json n> [output file "
                "(optional)]"
             << endl;
        return 0b00000001;
    }

    MPI_Init(&argc, &argv);
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    int aux = argc - 1;
    string out_file = argv[argc - 1];
    if (world_rank == 0) {
        if (argc == 2 || !ends_with(out_file, ".csv")) {
            const auto tp_utc{chrono::system_clock::now()};
            auto in_time_t = chrono::system_clock::to_time_t(tp_utc);
            stringstream ss;

            ss << put_time(localtime(&in_time_t), "%Y%m%d%H%M%S");
            out_file = "test_" + ss.str() + ".csv";
            cout << "Using default outfile: " << out_file << endl;
            aux++;
        }
    }
    const vector<string> fnames(argv + 1, argv + aux);
    //{
    //    "../saved_hypergraphs/hypergraphs_r/hypergraph_random_2025_04_28_19_39_28_6.json",
    //    "../saved_hypergraphs/hypergraphs_r/hypergraph_random_2025_04_28_19_39_28_7.json",
    //    "../saved_hypergraphs/hypergraphs_r/hypergraph_random_2025_04_28_19_39_28_8.json"
    //};
    run_multiple_test(out_file, fnames, world_rank, world_size);
    MPI_Finalize();

    return EXIT_SUCCESS;
}
