#pragma once
#include <iostream>
#include <vector>
#include <map>
#include <set>
#include <algorithm>
#include <numeric>
#include <random>
#include <cmath>
#include <string>
#include <sstream>
#include <stdexcept>
#include <limits>
#include <fstream>
#include <unordered_set>
#include <gmpxx.h>
#include "json.hpp"

using json = nlohmann::json;
using namespace std;

// Helper function templates
pair<vector<uint8_t>, int> int_to_bits(int x, int size);

set<int> vec_to_set(const vector<int>& v);

struct VectorIntCompare {
    bool operator()(const vector<int>& a, const vector<int>& b) const;
};

void print_vector_of_vectors(const vector<vector<int>>& F, string l);

template<typename T>
void print_vector(const vector<T>& v, string l) {
    cout << "Vector " << l ;
    int i =0;
    cout << "[";
    for (T val : v) {
        if (i++ > 0) cout << ", ";
        cout << val;
    }
    cout << "]" << endl;
}

string get_datetime_string();


class bipartite_graph {
public:
    map<int, vector<int>> Adj;
    map<int, vector<int>> Edges;
    set<int> Vertices;
    map<int, int> cnum;

    bipartite_graph(const vector<vector<int>>& H);
    bipartite_graph(const bipartite_graph& other);
    bipartite_graph& operator=(const bipartite_graph& other);

    vector<int> find_leaf_vertices();
    void add_edge(int u, int v);
    void add_edges(const vector<pair<int, int>>& E);
    vector<bipartite_graph> connected_components();
    void visit(int v_or_e, int current_cnum);
    bool intersection_property();
    void remove_set_vertex(const std::vector<int>& e);
    void merge_set_vertices() ;
    mpz_class count_hitting_set_brute_force(bool verbose = false);

    tuple<int, map<int, vector<int>>, map<int, vector<int>>, set<int>>
    compute_H_2(const vector<int>& p, int e_id);
    mpz_class compute_number_hitting_sets(mpz_class cut_off = 0);
    mpz_class count_hitting_set_recursive(mpz_class cut_off = 0);
    bool hitting_set(const vector<int>& s);
    pair<int, int> max_degree();
    int max_edge();
    size_t count_vertices();
    void print_graph(const string& label = "");
};

class hypergraph {
public:
    int n;
    set<vector<int>, VectorIntCompare> psi;
    vector<vector<uint8_t>> psi_values;
    vector<uint32_t> psi_dict;

    hypergraph(int num_vars = 0);

    pair<bool, mpz_class> check_dimension();
    void compute_values();
    long long self_dual();
    mpz_class self_dual_algorithms();
    void count_covered();
    void generate2(int k, const vector<int>& init_edge);
    void generate_all(int k);
    void generate_uniform();
    void generate_random(const vector<int>& init, int n_1, int n_2, int l);
    void generate_random_comb(int d_1, int d_2, int l);
    void generate_random_comb_new(int d_1, int d_2, long long l);
    void generate_random_nonuniform(int d_1, int d_2, long long l);
    void generate_random_exponential(int d_1, int d_2, long long l, double lambda);
    bool intersection(const vector<int>& s);
    vector<long long> get_number_list();
    int evaluate_f(const vector<uint8_t>& x_vec, int & w);
    mpz_class sum_f();
    pair<bool, vector<uint8_t>> search_x();
    bool check_counter_example(const vector<uint8_t>& x);
    size_t count_vertices();
    vector<vector<int>> get_ordered_list_psi();
    vector<int> min_length_hyperedge();
    vector<int> max_length_hyperedge();
    void save(const string& fname);
    void load(const string& fname);
    void load_uno_dat(const string& filename);
    map<string, double> stat();
    void print_psi();
    void print_bin_psi(bool compact = true);
};

void reduce_minimal(vector<vector<int>>& E);
bool algorithm_A(vector<vector<int>> F, vector<vector<int>> G);

// Template implementations must be in the header

template <typename T>
class CombinationGenerator {
public:
    CombinationGenerator(const vector<T>& pool, int r)
        : pool_(pool), r_(r), n_(pool.size()), done_(false)
    {
        combination_.resize(r_);
        if (r_ < 0 || r_ > n_) {
            done_ = true;
            return;
        }
        if (n_ == 0 && r_ == 0) {
            // Special case: one empty combination
            done_ = false;
            indices_.clear();
        } else if (n_ == 0 ) {
            done_ = true;
        } else {
            indices_.resize(r_);
            iota(indices_.begin(), indices_.end(), 0);
        }
        //build_combination();
        for (int idx =0; idx<r_; idx++) {
            combination_[idx]= pool_[indices_[idx]];
        }
        n_r = n_-r_;
    }

    // Returns  the current combination if not done, else empty vector
    const vector<T>& next() {
        if (done_) {
            result.clear();
            return result;
        }

        // Save current combination to return
        result = combination_;
        
        // Prepare indices for the next call
        int i;
        for (i = r_ - 1; i >= 0; --i) {
            if (indices_[i] != i +  n_r) {        //   (n_ - r_)) {
                break;
            }
            //combination_[i]= pool_[indices_[i]];
        }

        if (i < 0) {
            done_ = true;
        } else {
            indices_[i] += 1;
            combination_[i]= pool_[indices_[i]];
            for (int j = i + 1; j < r_; ++j) {
                indices_[j] = indices_[j - 1] + 1;
                combination_[j]= pool_[indices_[j]];
            }
            /*
            for (i--; i>=0; i--){
                combination_[i]= pool_[indices_[i]];
            }
            
            for (int idx =0; idx<r_; idx++) {
                combination_[idx]= pool_[indices_[idx]];
            }
            */
        }

        return result;
    }
    const vector<T>& get_current() const {
        return combination_;
    }

    void print_current() const {
        cout << "*** ";
        for (const T& elem : combination_) {
            cout << elem << " ";
        }
        cout << endl;
    }

    bool is_done() const {
        return done_;
    }

protected:
    void build_combination_old() {
        combination_.clear();
        for (int idx : indices_) {
            combination_.push_back(pool_[idx]);
        }
    }
    void build_combination() {
        
        for (int idx =0; idx<r_; idx++) {
            
            combination_[idx]= pool_[indices_[idx]];
            
        }
        //print_current();
    }

    const vector<T>& pool_;
    int r_;
    int n_;
    int n_r;
    vector<int> indices_;
    vector<T> combination_;
    vector<T> result;
    bool done_;
};


