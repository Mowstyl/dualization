#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <filesystem>
#include <chrono>
#include <iomanip>
#include <map>
#include <random>
#include <sstream>
#include "bipartite_hypergraph.hpp"
#include "bipartite_02.hpp"

using namespace std;

// Write a row to a CSV file
void write_csv_row(ofstream& file, const vector<string>& row) {
    for (size_t i = 0; i < row.size(); ++i) {
        file << row[i];
        if (i + 1 < row.size()) file << ";";
    }
    file << "\n";
}

// Implementation of the algorithm execution
vector<string> execute_algorithm(int algorithm_num, bipartite_graph& b, hypergraph& g) {
    string name;
    double elapsed = 0.0;
    string output = "";
    // Prepare data for algorithm_A: convert psi (set of vector) to vector of vector
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
        {
            auto search_result = g.search_x();
            ostringstream oss;
            oss << (search_result.first ? "self-dual" : "not self-dual");
            if (!search_result.first) {
                oss << " [";
                for (size_t i = 0; i < search_result.second.size(); ++i) {
                    oss << to_string(search_result.second[i]);
                    if (i + 1 < search_result.second.size()) oss << ",";
                }
                oss << "]";
            }
            output = oss.str();
        }
    } else if (algorithm_num == 3) {
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
                if (i>0) oss << ",";
                i++;
                oss << to_string(v);
            }
            oss << "]";
        }
        output = oss.str();
    } else if (algorithm_num ==4 ) {
        name = "sum_f";
        output = g.sum_f().get_str();
    }
    
    else {
        cerr << "Unknown algorithm number: " << algorithm_num << endl;
        return {};
    }
    
    auto end = chrono::high_resolution_clock::now();
    elapsed = chrono::duration<double>(end - start).count();

    cout << name << ": " << elapsed << " " << output << endl;
    return {name, to_string(elapsed), output};
}


void run_multiple_test(const string& fname, const vector<string>& hypergraphs) {

    ofstream csvFile(fname, ios::app);
    if (!csvFile.is_open()) {
        cerr << "Impossibile to open file " << fname << endl;
        return;
    }

    hypergraph g(0);
    for (const string& filepath : hypergraphs) {

        string filename = filesystem::path(filepath).filename().string();
        string datetime_string = get_datetime_string();
        g.load(filepath);
        int n = g.n;
        auto st = g.stat();
        bipartite_graph b(g.get_ordered_list_psi());
        cout << "hypergraph : " << filename << endl;
        cout << "Number of vertices: " << n << "; Number of hyperedges : " << b.Edges.size() << endl;

        vector<string> row_head = {
            filename,
            datetime_string,
            to_string(static_cast<int>(st["n"])),
            to_string(static_cast<int>(st["m"])),
            to_string(static_cast<int>(st["min_edge"])),
            to_string(static_cast<int>(st["max_edge"])),
            to_string(st["avg"])
        };
        
        for (int algorithm_num = 0; algorithm_num < 5; ++algorithm_num) {
            auto row = execute_algorithm(algorithm_num, b, g);
            vector<string> new_row = row_head;
            new_row.insert(new_row.end(), row.begin(), row.end());
            write_csv_row(csvFile, new_row);
        }

        cout << "----------------------------------------" << endl;

    }
    csvFile.close();
}



int main() {
    const vector<string> fname = {
        "../saved_hypergraphs/hypergraph_r/hypergraph_random_2025_04_28_19_39_28_6.json",
        "../saved_hypergraphs/hypergraph_r/hypergraph_random_2025_04_28_19_39_28_7.json",
        "../saved_hypergraphs/hypergraph_r/hypergraph_random_2025_04_28_19_39_28_8.json"
    };
    run_multiple_test("test.csv", fname);


    return 0;   
}