#include "bipartite_hypergraph.hpp"
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <string>
#include <thread>
#include <vector>

const string DIR_PATH = "./saved_hypergraphs/testing";

vector<int> generate_random_hyperedges(int n, int size) {
    /*
        Generates a random hyperedge with equal to size
        Returns a vector of integers representing the hyperedge.
    */
    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<> dis(0, n - 1);

    vector<uint8_t> e_bits(n, 0);
    fill_n(e_bits.begin(), size, 1); // Set the first j elements to 1
    shuffle(e_bits.begin(), e_bits.end(), gen);

    vector<int> s_vec;
    for (int bit_idx = 0; bit_idx < n; ++bit_idx) {
        if (e_bits[bit_idx] == 1) {
            s_vec.push_back(bit_idx);
        }
    }
    return s_vec;
}

pair<bool, mpz_class> check_dimension(hypergraph &g) {
    /*
        This function chek if condition (2.19) of Khachiyan paper is satisfied.
        It returns a pair containing a boolean indicating if the  sum
        $E_sum = \sum_{e \in \psi} 2^{n-|e|} is less than 2^{n-1}.
        and the sum of the sizes of the hyperedges and E_sum itself.
        PARAMETERS
        ------------------------
        g : hypergraph
            The hypergraph to check.
        RETURNS
        ------------------------
        A pair containing a boolean indicating if the condition is satisfied
        and the sum E_sum.
    */
    auto F = g.get_ordered_list_psi();
    mpz_class E_sum = 0;
    mpz_class M, E;
    mpz_ui_pow_ui(M.get_mpz_t(), 2, g.n - 1); // M = 2^(n-1)
    for (const auto &e : F) {
        mpz_ui_pow_ui(E.get_mpz_t(), 2, g.n - e.size()); // E = 2^(n - |e|)
        E_sum += E;
    }
    // E_sum = \sum_{e \in F} 2^{n-|e|}
    return {E_sum < M, E_sum};
}
void generate_random_save_comb(int V_min, int V_max, int step, int range,
                               string prefix) {
    /* *
    This function generates and saves random hypergraphs using the
    generate_random_comb method. For each n in [V_min, V_max] with the given
    step, it creates a hypergraph with n vertices. The number of iterations is
    2^(n+2) (exponential growth). Random hyperedges of sizes in [range, n-range]
    are generated and the hypergraph is saved as a JSON file named with a
    timestamp and n in DIR_PATH. After each generation, the function prints the
    elapsed time and pauses for 1 second. PARAMETERS
    ------------------------
    V_min : int
        Minimum number of vertices in the hypergraph.
    V_max : int
        Maximum number of vertices in the hypergraph.
    step : int
        Step size for the number of vertices.
    range : int
        Range for the size of the hyperedges.
    prefix : string
        Prefix for the filename of the saved hypergraph.
    * */
    int n_iter_exp = 2;

    long long n_iter;

    for (int n = V_min; n <= V_max; n += step) {
        // Get current datetime as string
        string datetime_string = get_datetime_string();
        cout << "Generating hypergraph with " << n << " vertices..." << endl;

        n_iter = 1L << (n + n_iter_exp); // Default value for n_iter

        hypergraph g(n);
        auto ini = chrono::high_resolution_clock::now();
        g.generate_random_comb_new(range, n - range, n_iter);
        auto fin = chrono::high_resolution_clock::now();
        double elapsed_time = chrono::duration<double>(fin - ini).count();

        ostringstream fname_oss;
        fname_oss << prefix << datetime_string << "_" << n << ".json";
        string fname = DIR_PATH + "/" + fname_oss.str();

        g.save(fname);
        cout << "Time: " << elapsed_time << endl;
        cout << "Hyperedges: " << g.psi.size() << endl;
        cout << "---------------------" << endl;
        this_thread::sleep_for(chrono::seconds(1));
    }
}

void generate_random_save_comb01(int V_min, int V_max, int step,
                                 string prefix) {
    /*
    This function generates and saves random hypergraphs using the
    generate_random_exponential method. For each n in [V_min, V_max] with the
    given step, it creates a hypergraph with n vertices. The number of
    iterations is  n^2 (polynomial growth). The hypergraph is initialized with a
    hyperedge of size m, where m is the maximum size of the minimum size
    hyperedge. he value $m$ is determined as the logarithm base 2 of n, rounded
    to the nearest integer. Random hyperedges of sizes in [m, n-m] are generated
    using an exponential distribution with a rate parameter (λ) of 2. If the
    condition (2.19) of the Khachiyan paper is not satisfied, another set of
    hyperedges is generated but with a different value of λ, which is increased
    by 0.1 in each iteration. The process continues until the condition is
    satisfied. After the initial hypergraph is generated,  another 4*n^2
    iterations are executed inorder to add another set of hyperedges. this time
    the parameter λ is set to 2.0. After each generation, the function prints
    the elapsed time and pauses for 1 second.

    PARAMETERS
    ------------------------
    V_min : int
        Minimum number of vertices in the hypergraph.
    V_max : int
        Maximum number of vertices in the hypergraph.
    step : int
        Step size for the number of vertices.
    prefix : string
        Prefix for the filename of the saved hypergraph.

    */
    long long n_iter;
    double lambda = 3.0;

    for (int n = V_min; n <= V_max; n += step) {
        // Get current datetime as string
        string datetime_string = get_datetime_string();
        cout << "Generating hypergraph with " << n << " vertices... ";

        n_iter = n * n; // We use n^2 for polynomial growth

        hypergraph g(n);
        int m = static_cast<int>(
            round(log2(n))); // Maximum size of minimum size hyperedge

        vector<int> e_0 = generate_random_hyperedges(
            n, m); // Initial hyperedge with logaritmic size

        auto ini_all = chrono::high_resolution_clock::now();
        int iter = 0;
        int n_try = 0;
        auto r = g.check_dimension();
        while (r.first) {

            g.psi.clear();
            g.psi.insert(e_0); // Insert initial hyperedge
            g.generate_random_exponential(e_0.size(), n - e_0.size(), n_iter,
                                          lambda);

            iter++;
            if (iter % 50 == 0) {

                cout << "Iteration: " << iter;

                lambda += 0.1; // Increase the delta  parameter
                cout << "e_0 size: " << e_0.size() << " lambda: " << lambda
                     << endl;
            }

            r = g.check_dimension();
        }
        cout << "Hyperedges before the addition: " << g.psi.size() << endl;
        // continue to generate another 4 rounds of other hyperedges
        // but with lambda = 2
        double lambda1 = 2.0;
        for (int i = 0; i < 4; ++i) {
            g.generate_random_exponential(e_0.size(), n - e_0.size(), n_iter,
                                          lambda1);
        }
        cout << "Hyperedges after the 2cond addition: " << g.psi.size() << endl;
        auto fin_all = chrono::high_resolution_clock::now();
        double elapsed_time =
            chrono::duration<double>(fin_all - ini_all).count();

        ostringstream fname_oss;
        fname_oss << prefix << datetime_string << "_" << n << ".json";
        string fname = DIR_PATH + "/" + fname_oss.str();

        g.save(fname);
        cout << "Time: " << elapsed_time << endl;
        cout << "Hyperedges final: " << g.psi.size() << endl;
        this_thread::sleep_for(chrono::seconds(1));
    }
}

int main(int argc, char *argv[]) {
    if (argc != 5) {
        cout << "Syntax: " << argv[0] << " <V_min> <V_max> <step> <prefix>"
             << endl;
        return 0b00000001;
    }

    // generate_random_save_comb01(210, 260, 10, "hypergraph_random_poly_4n2_");
    generate_random_save_comb01(atoi(argv[1]), atoi(argv[2]), atoi(argv[3]),
                                argv[4]);

    return 0;
}