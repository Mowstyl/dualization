#include "bipartite_hypergraph.hpp"
#include <omp.h>

pair<vector<uint8_t>, int> mpz_to_bits(const mpz_class &num, int size) {

    /*
        Converts an mpz_class number to a vector of bits of size `size`.
        Returns a pair containing the vector of bits and the number of 1s in it.
    */
    int w = 0;
    vector<uint8_t> bits(size, 0);
    mpz_class n = num;
    for (int i = 0; i < size; ++i) {
        bits[i] = mpz_tstbit(n.get_mpz_t(), i);
        if (bits[i] == 1) {
            w++;
        }
    }
    return {bits, w};
}

pair<vector<uint8_t>, int> int_to_bits(int x, int size) {
    /*
        Converts an integer x into a vector of bits of size `size`.
        Returns a pair containing the vector of bits and the number of 1s in it.
    */
    vector<uint8_t> bits(size);
    int pos = 0;
    int w = 0;
    while (x > 0) {
        bits[pos] = (x & 1);
        if (bits[pos] == 1) {
            w++;
        }
        x = x >> 1;
        pos++;
    }
    return {bits, w};
}

// Implementations for VectorIntCompare and vec_to_set
bool VectorIntCompare::operator()(const vector<int> &a,
                                  const vector<int> &b) const {
    if (a.size() != b.size()) {
        return a.size() < b.size();
    }
    return a < b;
}

set<int> vec_to_set(const vector<int> &v) {
    return set<int>(v.begin(), v.end());
}

bipartite_graph::bipartite_graph(const vector<vector<int>> &H) {
    if (H.empty()) {
        return;
    }

    // Populate Vertices set
    for (const auto &e : H) {
        for (int v : e) {
            Vertices.insert(v);
        }
    }

    // Find max vertex ID to assign edge IDs
    int max_v = 0;
    if (!Vertices.empty()) {
        max_v = *max_element(Vertices.begin(), Vertices.end());
    }
    max_v++;

    // Populate Edges dictionary
    for (size_t i = 0; i < H.size(); ++i) {
        Edges[i + max_v] = H[i];
    }

    // Build adjacency lists for the bipartite graph
    for (const auto &pair_e : Edges) {
        int edge_id = pair_e.first;
        const vector<int> &vertices_in_edge = pair_e.second;
        for (int v : vertices_in_edge) {
            add_edge(v, edge_id);
        }
    }
}

bipartite_graph::bipartite_graph(const bipartite_graph &other)
    : Adj(other.Adj), Edges(other.Edges), Vertices(other.Vertices) {}

bipartite_graph &bipartite_graph::operator=(const bipartite_graph &other) {
    if (this != &other) {
        Adj = other.Adj;
        Edges = other.Edges;
        Vertices = other.Vertices;
    }
    return *this;
}

void bipartite_graph::add_edge(int u, int v) {
    if (Adj.find(v) == Adj.end()) {
        Adj[v] = {u};
    } else {
        Adj[v].push_back(u);
    }
    if (Adj.find(u) == Adj.end()) {
        Adj[u] = {v};
    } else {
        Adj[u].push_back(v);
    }
}

void bipartite_graph::add_edges(const vector<pair<int, int>> &E) {
    for (const auto &edge : E) {
        add_edge(edge.first, edge.second);
    }
}

vector<int> bipartite_graph::find_leaf_vertices() {
    vector<int> leaves;
    for (const auto &v : Vertices) {
        if (Adj[v].size() == 1) {
            leaves.push_back(v);
        }
    }
    return leaves;
}

vector<bipartite_graph> bipartite_graph::connected_components() {
    cnum.clear();
    for (int v : Vertices) {
        cnum[v] = 0;
    }
    for (const auto &pair_e : Edges) {
        cnum[pair_e.first] = 0;
    }

    vector<bipartite_graph> CC;
    int current_cnum = 1;

    for (int v : Vertices) {
        if (cnum[v] == 0) {
            visit(v, current_cnum);
            bipartite_graph C({}); // Create an empty graph for the component
            for (const auto &adj_pair : Adj) {
                if (cnum[adj_pair.first] == current_cnum) {
                    C.Adj[adj_pair.first] =
                        adj_pair.second; // Copy adjacency list

                    if (Vertices.count(adj_pair.first)) { // If it's a vertex
                        C.Vertices.insert(adj_pair.first);
                    } else { // If it's an edge
                        C.Edges[adj_pair.first] = Edges[adj_pair.first];
                    }
                }
            }
            CC.push_back(C);
            current_cnum++;
        }
    }
    cnum.clear(); // Clear temporary cnum
    return CC;
}
void bipartite_graph::visit(int v_or_e, int current_cnum) {
    cnum[v_or_e] = current_cnum;
    // Check if Adj[v_or_e] exists before iterating
    if (Adj.count(v_or_e)) {
        for (int neighbor : Adj[v_or_e]) {
            if (cnum[neighbor] == 0) {
                visit(neighbor, current_cnum);
            }
        }
    }
}

bool bipartite_graph::intersection_property() {
    for (const auto &pair_e : Edges) {
        map<int, bool> tested;
        for (const auto &other_pair_e : Edges) {
            tested[other_pair_e.first] = false;
        }
        tested[pair_e.first] = true;

        int edge_id = pair_e.first;
        if (Adj.count(edge_id)) {
            for (int v : Adj[edge_id]) {
                if (Adj.count(v)) {
                    for (int f : Adj[v]) {
                        tested[f] = true;
                    }
                }
            }
        }
        for (const auto &test_pair : tested) {
            if (!test_pair.second) {
                return false;
            }
        }
    }
    return true;
}

void bipartite_graph::remove_set_vertex(const vector<int> &e) {
    /*
        Input: vector<int> e
        Removes a set vertex e from the bipartite graph.
        For each vertex v adjacent to e, removes v from all edges containing it,
        removes v from the adjacency lists, and from the vertex set.
        Finally, removes the adjacency and edge for e itself.
    */

    bool found = false;
    int ke = -1;
    // Find the key ke such that Edges[ke] == e
    for (const auto &pair : Edges) {
        if (pair.second == e) {
            found = true;
            ke = pair.first;
            break;
        }
    }
    if (!found) {
        cout << "node {";
        for (size_t i = 0; i < e.size(); ++i) {
            cout << e[i];
            if (i + 1 < e.size())
                cout << ",";
        }
        cout << "} not present in the graph" << endl;
        return;
    }

    // Copy the adjacency list of ke (vertices adjacent to e)
    vector<int> v_remove = Adj[ke];

    for (int v : v_remove) {
        // For each edge ie adjacent to v, remove v from Edges[ie] and from
        // Adj[ie]
        vector<int> adj_edges = Adj[v]; // Copy to avoid iterator invalidation
        for (int ie : adj_edges) {
            // Remove v from Edges[ie]
            auto &edge_vec = Edges[ie];
            edge_vec.erase(remove(edge_vec.begin(), edge_vec.end(), v),
                           edge_vec.end());
            // Remove v from Adj[ie]
            auto &adj_vec = Adj[ie];
            adj_vec.erase(remove(adj_vec.begin(), adj_vec.end(), v),
                          adj_vec.end());
        }
        Adj.erase(v);      // Remove adjacency of v
        Vertices.erase(v); // Remove v from vertex set
    }
    Adj.erase(ke);   // Remove adjacency of e
    Edges.erase(ke); // Remove the edge e
}

void bipartite_graph::merge_set_vertices() {
    /*
        After deleting a set vertex from the bipartite graph, the hypergraph can
       have multiple edges. This function merges such duplicate edges by
       removing all but one instance of each duplicate edge.

        Returns
        -------
        None.
    */

    // Collect all edge keys
    vector<int> keys;
    for (const auto &pair : Edges) {
        keys.push_back(pair.first);
    }

    // Find duplicate edges and store their keys in edges_removed
    set<int> edges_removed;
    for (size_t i = 0; i + 1 < keys.size(); ++i) {
        for (size_t j = i + 1; j < keys.size(); ++j) {
            if (Edges[keys[i]] == Edges[keys[j]]) {
                edges_removed.insert(keys[j]);
            }
        }
    }

    // Remove duplicate edges from adjacency lists and from Edges/Adj
    for (int e : edges_removed) {
        for (int v : Edges[e]) {
            auto &adj_vec = Adj[v];
            adj_vec.erase(remove(adj_vec.begin(), adj_vec.end(), e),
                          adj_vec.end());
        }
        Adj.erase(e);
        Edges.erase(e);
    }
}

mpz_class bipartite_graph::count_hitting_set_brute_force(bool verbose) {
    mpz_class n_hitting_sets = 0;
    vector<int> vertices_vec(Vertices.begin(), Vertices.end());
    vector<int> s;

    for (size_t i = 1; i <= vertices_vec.size(); ++i) {
        for (CombinationGenerator<int> gen(vertices_vec, i); !gen.is_done();) {
            s = gen.next();
            if (hitting_set(s)) {
                n_hitting_sets++;
            }
        }
    }
    return n_hitting_sets;
}

tuple<int, map<int, vector<int>>, map<int, vector<int>>, set<int>>
bipartite_graph::compute_H_2(const vector<int> &p, int e_id) {
    /*
            Parameters
        ----------
        p : list or tuple
            the set of vertices neede to be removed .
        e : int
            index of the edges considered.
        Computes $H_2$ as reported in the paper. First $H_1=H \setmins p$
        is computed. Then
        Returns
        -------
        n_i: int
            the set of vertices in $V(H_1)\setminus (V(H_2) \cup e\setminus p)$
        edges_backup:
        adj_backup:
        vertices_backup:
            the backup of the edges, adjacency lists and vertices before the
       deletion.
    */
    // Make a copy of the edge e (which is Edges[e_id])
    vector<int> E_copy = Edges[e_id];

    map<int, vector<int>> edges_backup;
    map<int, vector<int>> adj_backup;
    set<int> vertices_backup;

    // Backup current vertices in E_copy for later restoration
    for (int v : E_copy) {
        vertices_backup.insert(v);
    }

    /* Compute H_1 = H \setminus p */
    for (int v : p) {
        if (Adj.count(v)) {
            adj_backup[v] = Adj[v]; // Backup the adj list of v (if not already
                                    // backed up by p)
            for (int current_e_id : Adj[v]) {
                // Backup edges modified and their adjacency lists
                if (Edges.count(current_e_id)) {
                    if (edges_backup.find(current_e_id) == edges_backup.end()) {
                        edges_backup[current_e_id] = Edges[current_e_id];
                        adj_backup[current_e_id] = Adj[current_e_id];
                    }
                    // Remove v from the edge
                    auto &edge_vec = Edges[current_e_id];
                    // the line below wa replaced with the commented one
                    edge_vec.erase(find(edge_vec.begin(), edge_vec.end(), v));
                    // edge_vec.erase(remove(edge_vec.begin(), edge_vec.end(),
                    // v), edge_vec.end());

                    // Remove v from the adjacency list of the edge
                    auto &adj_vec = Adj[current_e_id];
                    // the line below wa replaced with the commented one
                    adj_vec.erase(find(adj_vec.begin(), adj_vec.end(), v));
                    // adj_vec.erase(remove(adj_vec.begin(), adj_vec.end(), v),
                    // adj_vec.end());
                }
            }
        }
        // Delete v from adjacency lists and from Vertices
        Adj.erase(v);
        Vertices.erase(v);
        // Remove v from E_copy (since it's now s)
        E_copy.erase(find(E_copy.begin(), E_copy.end(), v));
        // E_copy.erase(remove(E_copy.begin(), E_copy.end(), v), E_copy.end());
    }

    /* Compute H_2 = H_1 \setminus N_{H_1}(s) */
    set<int> deleted_edges;
    for (int v : E_copy) { // E_copy now represents s
        if (Adj.count(v)) {
            adj_backup[v] = Adj[v]; // Backup the adj list of v
            for (int current_e_id :
                 Adj[v]) { // Backup adjacency list of edges to be deleted
                if (edges_backup.find(current_e_id) == edges_backup.end()) {
                    edges_backup[current_e_id] = Edges[current_e_id];
                    adj_backup[current_e_id] = Adj[current_e_id];
                }
                deleted_edges.insert(current_e_id);
            }
        }
    }

    for (int e_v_id : deleted_edges) { // For each deleted edge
        if (Adj.count(e_v_id)) {
            vector<int> adj_e_v_copy = Adj[e_v_id]; // Copy to iterate safely
            for (int u : adj_e_v_copy) { // For each vertex u in the adjacency
                                         // list of e_v_id
                if (adj_backup.find(u) ==
                    adj_backup.end()) { // Backup the vertex's u adjacency list
                    adj_backup[u] = Adj[u];
                    vertices_backup.insert(u);
                }
                if (Adj.count(u)) { // delete the edge e_v_id from the adjacency
                                    // list of u
                    auto &adj_vec = Adj[u];

                    adj_vec.erase(find(adj_vec.begin(), adj_vec.end(), e_v_id));
                    // adj_vec.erase(remove(adj_vec.begin(), adj_vec.end(),
                    // e_v_id), adj_vec.end());
                }
            }
            // Now delete the  adjacency list of e_v_id and the edge itself
            Adj.erase(e_v_id);
            Edges.erase(e_v_id);
        }
    }

    /* Delete the vertices in E_copy (s) */
    for (int v : E_copy) {
        Adj.erase(v);
        Vertices.erase(v);
    }

    /* Count the vertices in V(H_1) \setminus (V(H_2) \cup E_copy) and delete
     * them */
    vector<int> deleted_vertices_isolated;
    int n_i = 0;
    for (int v : Vertices) {
        if (!Adj.count(v) ||
            Adj[v].empty()) { // Check if adjacency list exists and is empty
            n_i++;
            deleted_vertices_isolated.push_back(v);
        }
    }

    for (int v : deleted_vertices_isolated) {
        Vertices.erase(v);
        Adj.erase(v);
    }

    return make_tuple(n_i, edges_backup, adj_backup, vertices_backup);
}

mpz_class bipartite_graph::compute_number_hitting_sets(mpz_class cut_off) {
    /*
        Computes the number of hitting sets recursively.
        If cut_off is set, it will stop recursion when the number of hitting
       sets exceeds cut_off.
    */
    mpz_class n_hit = 1;
    if (Edges.empty()) {
        return n_hit;
    }
    if (Edges.size() == 1) {
        // If there is only one edge, the hitting set is the vertices in that
        // edge
        mpz_ui_pow_ui(n_hit.get_mpz_t(), 2, Edges.begin()->second.size());
        return n_hit - 1;
    }

    vector<bipartite_graph> CC = connected_components();
    // mpz_class n_hit = 1;

    for (auto &BH_i : CC) {
        n_hit *= BH_i.count_hitting_set_recursive();
    }
    return n_hit;
}
mpz_class bipartite_graph::count_hitting_set_recursive(mpz_class cut_off) {
    mpz_class n_hit = 0;
    int min_len = numeric_limits<int>::max();
    int e_min_len_id = -1; // The id of the edge with minimum cardinality

    for (const auto &pair_e : Edges) {
        if (pair_e.second.size() < min_len) {
            e_min_len_id = pair_e.first;
            min_len = pair_e.second.size();
        }
    }

    // This if should be deleted
    if (e_min_len_id == -1) { // Should not happen if Edges is not empty
        return 0;
    }

    vector<int> E_current_edge = Edges[e_min_len_id];
    vector<int> elements_for_combinations(E_current_edge.begin(),
                                          E_current_edge.end());
    vector<int> p;
    for (int i = 0; i < min_len; ++i) {
        for (CombinationGenerator<int> gen(elements_for_combinations, i);
             !gen.is_done();) {
            p = gen.next();

            auto [n_i, e_b, a_b, v_b] = this->compute_H_2(p, e_min_len_id);

            mpz_class n_r = this->compute_number_hitting_sets(cut_off);

            mpz_class n_r1;
            mpz_ui_pow_ui(n_r1.get_mpz_t(), 2, n_i);
            n_hit += n_r1 * n_r;

            // Revert changes back to *this for the next iteration
            // This is done implicitly because we are working on temp_graph
            // and then restoring values to 'this' from the backups
            // However, the Python code applies changes directly to 'self', then
            // restores To mimic Python's behavior, we need to apply the backup
            // directly to 'this'
            for (const auto &pair : e_b) {
                this->Edges[pair.first] = pair.second;
            }
            for (const auto &pair : a_b) {
                this->Adj[pair.first] = pair.second;
            }
            for (int v : v_b) {
                this->Vertices.insert(v);
            }
        }
    }
    return n_hit;
}

bool bipartite_graph::hitting_set(const vector<int> &s) {
    unordered_set<int> hit;
    int count = 0;
    for (int v : s) {
        auto it = Adj.find(v);
        if (it != Adj.end()) {
            for (int e : it->second) {
                // Insert returns a pair; if insertion happened, increment count
                if (hit.insert(e).second) {
                    count++;
                    // Early exit if all edges are hit
                    if (count == Edges.size())
                        return true;
                }
            }
        }
    }
    return count == Edges.size();
}

pair<int, int> bipartite_graph::max_degree() {
    int max_deg = 0;
    int vertex_with_max_deg = 0;
    for (int v : Vertices) {
        if (Adj.count(v) && Adj[v].size() > max_deg) {
            max_deg = Adj[v].size();
            vertex_with_max_deg = v;
        }
    }
    return {vertex_with_max_deg, max_deg};
}

int bipartite_graph::max_edge() {
    int m = 0;
    for (const auto &pair_e : Edges) {
        if (pair_e.second.size() > m) {
            m = pair_e.second.size();
        }
    }
    return m;
}

void bipartite_graph::print_graph(const string &label) {
    cout << label;
    for (int v : Vertices) {
        if (Adj.count(v) && Adj[v].empty()) {
            cout << v << endl;
        } else if (Adj.count(v)) {
            for (int e : Adj[v]) {
                cout << v << " " << e << " {";
                if (Edges.count(e)) {
                    for (size_t i = 0; i < Edges[e].size(); ++i) {
                        cout << Edges[e][i];
                        if (i < Edges[e].size() - 1) {
                            cout << ",";
                        }
                    }
                }
                cout << "}" << endl;
            }
        }
    }
    for (const auto &pair_e : Edges) {
        cout << "edge " << pair_e.first << " {";
        if (Adj.count(pair_e.first)) {
            for (size_t i = 0; i < Adj[pair_e.first].size(); ++i) {
                cout << Adj[pair_e.first][i];
                if (i < Adj[pair_e.first].size() - 1) {
                    cout << ",";
                }
            }
        }
        cout << "}" << endl;
    }
}
// End of bipartite_graph implementation

// Implementation of hypergraph class
hypergraph::hypergraph(int num_vars) : n(num_vars) {}

pair<bool, mpz_class> hypergraph::check_dimension() {
    /*
        Checks if the condition (2.1) Fredman and  Khachiyan paper.
        Returns a pair containing a boolean indicating if the hypergraph is of
       dimension 2^n-1 and the sum of the sizes of the hyperedges.
    */
    auto F = this->get_ordered_list_psi();
    mpz_class E_sum = 0;
    mpz_class M, E;
    mpz_ui_pow_ui(M.get_mpz_t(), 2, n - 1);
    for (const auto &e : F) {
        mpz_ui_pow_ui(E.get_mpz_t(), 2, n - e.size());
        // E_sum += pow(2, -static_cast<double>(e.size()));
        E_sum += E;
    }

    return {E_sum < M, E_sum};
}

void hypergraph::generate_random_exponential(int d_1, int d_2, long long l,
                                             double lambda) {
    /**
     * Generates a random hypergraph with hyperedges whose sizes follow
     * a distribution derived from an exponential distribution, and
     * ensures the generated hyperedges satisfy an intersection property.
     *
     * This function populates the `psi` member (presumably a set of hyperedges,
     * where each hyperedge is a vector of vertex IDs) based on random
     * generation. The size of each hyperedge is determined by sampling from an
     * exponential distribution and then scaling/clamping the result.
     *
     * PARAMETRES
     * ------------------------
     * d_1 An integer representing the minimum desired size for a hyperedge,
     * or a lower bound for the scaled exponential value.
     *
     * d_2 An integer representing the maximum desired size for a hyperedge,
     * or an upper bound for the scaled exponential value.
     *
     * l A long long integer representing the number of iterations or
     * attempts to generate hyperedges. The final number of hyperedges
     * in `psi` might be less than `l` due to the intersection check.
     *  lambda A double representing the rate parameter (λ) for the
     * exponential distribution. A higher lambda means values are
     * more concentrated near 0.
     *
     * OUTPUT
     * ------------------------
     * This function does not return a value. Instead, it modifies the
     * `hypergraph` object it is called on, specifically the `psi` member.
     * The `psi` member of the `hypergraph` object is modified. It will contain
     * a set of unique hyperedges that satisfy the `intersection` property.
     * Each hyperedge is a `std::vector<int>` representing the vertices it
     * contains. The `n` member of the `hypergraph` object (number of vertices)
     * is used internally for generating hyperedges of appropriate size and
     * creating bit vectors.
     */
    double value;
    int j;
    static thread_local mt19937 *generator = nullptr;
    if (!generator) {
        generator =
            new mt19937(chrono::system_clock::now().time_since_epoch().count() +
                        omp_get_thread_num());
    }
    exponential_distribution<> distribution(lambda);

    for (long long i = 0; i < l; ++i) {
        value = distribution(*generator); // Generate a random value from the
                                          // exponential distribution
        j = static_cast<int>(round(value / 5 * (n - 2 * d_1))) +
            d_1; // Round to the nearest integer
        if (j > d_2) {
            j = d_2; // Ensure j does not exceed d_2
        }
        // cout << "j: " << j << endl;
        vector<uint8_t> e_bits(n, 0);
        fill_n(e_bits.begin(), j, 1); // Set the first j elements to 1
        shuffle(e_bits.begin(), e_bits.end(), *generator);

        vector<int> s_vec;
        for (int bit_idx = 0; bit_idx < n; ++bit_idx) {
            if (e_bits[bit_idx] == 1) {
                s_vec.push_back(bit_idx);
            }
        }

        // sort(s_vec.begin(), s_vec.end()); // Ensure the vector is sorted

        if (intersection(s_vec)) {
            psi.insert(s_vec);
        }
    }
}

bool hypergraph::intersection(const vector<int> &s) {
    /*
        Determine if the intersection property is satisfied for s with
        respect to this->psi. That is, if s intersects every set in this->psi
       and does not contains any set in this->psi and is not contained in any
        set in this->psi.
        The set this->psi must be sorted by size and lexicographically .

        Parameters
        ----------
        s : iterable
            the set to be tested. It must be sorted in ascending order.

        Returns
        -------
        bool
            True if the intersection property is satisfied. False otherwise
    */
    size_t ls = s.size();

    for (const auto &t : psi) {
        // Both s and t are sorted
        size_t i = 0, j = 0, count = 0;
        while (i < ls && j < t.size()) {
            if (s[i] < t[j]) {
                ++i;
            } else if (s[i] > t[j]) {
                ++j;
            } else {
                ++count;
                ++i;
                ++j;
            }
        }
        if (count == 0 || count == ls || count == t.size()) {
            return false;
        }
    }
    return true;
}

mpz_class hypergraph::sum_f() {
    /*
        This function computes the sum of 1-f(x) for all x in {0,1}^n. That is
        it computes the number of vectors x such that f(x) = 0.
        It uses the evaluate_f function to determine if f(x) = 0 or 1.

        It is mandatory that psi is ordered by length and lexicographically
        for this function to properly work.
    */
    mpz_class s = 0;
    vector<uint8_t> x_vec;
    mpz_class w;
    mpz_ui_pow_ui(w.get_mpz_t(), 2, n); // w = 2^n
    for (mpz_class x = 0; x < w; x++) {
        auto res = mpz_to_bits(x, n);
        x_vec = res.first;  // Get the bits as a vector
        int w = res.second; // Get the weight of the vector
        if (evaluate_f(x_vec, w) == 0) {
            s++;
        }
    }
    return s; // Return the number of vectors x such that f(x) = 0
}

int hypergraph::evaluate_f(const vector<uint8_t> &x_vec, int &w) {
    /*
        This function evaluates the function f(x) for a given binary vector
       x_vec. It checks if there is at least one hyperedges in psi that is
       satisfied by the vector x_vec. If one term is satisfied, it returns 1,
       otherwise it returns 0. It is mandatory that psi is ordered by length
       and,  for hyperedges of equal length, they are ordered lexicographically
        for this function to properly work.

        Parameters
        ----------
        x_vec : list of uint8_t
            The binary vector representing the input x.
        w : int
            The hamming weight of x_vec, i.e., the number of 1s in x_vec.

        Returns
        -------
        int
            Returns 1 if at least one hyperedges is satisfied, otherwise returns
       0.
    */
    for (const vector<int> &e : psi) {
        /*Since psi is ordered by length, if e.size() > w, we can return 0
        This is because if the hyperedge e has more vertices than the number of
        1s in x_vec, it cannot be satisfied by x_vec. For example, if e = {x_1,
        x_3, x_4, x_5} and x_vec = [0, 1, 0, 1, 1, 0], then w = 3 and e.size() =
        4 In this case, we can immediately return 0 because x_vec cannot satisfy
        e. Since psi is ordered by length, we can stop checking further
        hyperedges as they will also have size greater than w.

        */
        if (e.size() > w)
            return 0;
        int s = 0;
        for (int v : e) {
            if (x_vec[v] ==
                1) { // e = {x_1, x_3, x_4, x_5} and x_vec = [0, 1, 0, 1, 1, 0]
                s++;
            }
        }
        if (s == e.size()) {
            /*
            cout << "e: "  ;
            for (int v : e) {
                cout << v << " ";
            }
            cout << "; x_vec: ";
            for (int i = 0; i < x_vec.size(); ++i) {
                cout << (int)x_vec[i];
            }
            cout << "; w: " << w << endl;
            */
            return 1;
        }
    }

    return 0;
}

pair<bool, vector<uint8_t>> hypergraph::search_x() {
    /* This function searches for a binary vector x of length n such that
        f(x) = f(\overline{x}) = 0, where \overline{x} is the complement of x.
        It returns a pair containing a boolean indicating if such x exists and
       the vector x itself. If no such x is found, it returns true and an empty
       vector.

        The function iterates through all combinations of vertices in the
       hypergraph and checks if the intersection property holds. If it finds a
       valid x, it returns it.

        Parameters
        ----------
        None

        Returns
        -------
        pair<bool, vector<uint8_t>>
            A pair where the first element is true if no such x is found,
            and the second element is the binary vector x if found.
            If no such x is found, the second element is an empty vector.
        The function assumes that psi is ordered by length and
       lexicographically. It is mandatory that psi is ordered by length and
       lexicographically for this function to properly work. The function uses
       the evaluate_f function to check if f(x) = 0 or 1. It also uses the
       min_length_hyperedge function to determine the minimum length of
        hyperedges in psi, which is used to optimize the search.

    */

    int min_len = 0;
    if (!psi.empty()) {
        min_len = min_length_hyperedge().size();
    }

    int w = n / 2;

    vector<int> range_n(n);
    iota(range_n.begin(), range_n.end(), 0);
    vector<int> v_indices;

    for (int i = 1; i <= w; ++i) {
        for (CombinationGenerator<int> gen(range_n, i); !gen.is_done();) {
            v_indices = gen.next();

            // Check if n is even and if we are at the half-way point.
            if (n % 2 == 0 && i == w && v_indices[0] != 0) {
                break;
            }

            vector<uint8_t> x(n, 0);
            for (int idx : v_indices) {
                x[idx] = 1;
            }

            if (i < min_len) {
                vector<uint8_t> one_minus_x(n);
                for (int k = 0; k < n; ++k) {
                    one_minus_x[k] = 1 - x[k];
                }
                int w_val = n - i;
                if (evaluate_f(one_minus_x, w_val) == 0) {
                    // if (evaluate_f_old(one_minus_x) == 0) {
                    return {false, x};
                }
            } else {
                if (evaluate_f(x, i) == 0) {
                    // if (evaluate_f_old(x) == 0) {
                    vector<uint8_t> one_minus_x(n);
                    for (int k = 0; k < n; ++k) {
                        one_minus_x[k] = 1 - x[k];
                    }
                    int w_val = n - i;
                    if (evaluate_f(one_minus_x, w_val) == 0) {
                        // if (evaluate_f_old(one_minus_x) == 0) {
                        return {false, x};
                    }
                }
            }
        }
    }
    return {true, {}}; // Return true and an empty vector if no such x is found
}

unsigned long long combinations(int n, int k) {
    if (k > n)
        return 0;
    if (k * 2 > n)
        k = n - k;
    if (k == 0)
        return 1;

    unsigned long long result = n;
    for (int i = 2; i <= k; ++i) {
        result *= (n - i + 1);
        result /= i;
    }

    return result;
}

vector<int> ith_combination(int n, int k, unsigned long long i) {
    vector<int> result;
    int current = 0;
    unsigned long long count;

    for (int remaining = k; remaining > 0; ++current) {
        count = combinations(n - 1 - current, remaining - 1);

        if (i < count) {
            result.push_back(current);
            --remaining;
        } else {
            i -= count;
        }
    }
    assert((int)result.size() == k);
    return result;
}

pair<bool, vector<uint8_t>> hypergraph::search_x_par() {
    /* This function searches for a binary vector x of length n such that
        f(x) = f(\overline{x}) = 0, where \overline{x} is the complement of x.
        It returns a pair containing a boolean indicating if such x exists and
       the vector x itself. If no such x is found, it returns true and an empty
       vector.

        The function iterates through all combinations of vertices in the
       hypergraph and checks if the intersection property holds. If it finds a
       valid x, it returns it.

        Parameters
        ----------
        None

        Returns
        -------
        pair<bool, vector<uint8_t>>
            A pair where the first element is true if no such x is found,
            and the second element is the binary vector x if found.
            If no such x is found, the second element is an empty vector.
        The function assumes that psi is ordered by length and
       lexicographically. It is mandatory that psi is ordered by length and
       lexicographically for this function to properly work. The function uses
       the evaluate_f function to check if f(x) = 0 or 1. It also uses the
       min_length_hyperedge function to determine the minimum length of
        hyperedges in psi, which is used to optimize the search.

    */

    int min_len = 0;
    if (!psi.empty()) {
        min_len = min_length_hyperedge().size();
    }

    int w = n / 2;

    vector<int> range_n(n);
    iota(range_n.begin(), range_n.end(), 0);

    for (int i = 1; i <= w; ++i) {
        volatile bool flag = false;
        volatile bool exit = false;
        vector<uint8_t> out_x(n, 0);
        unsigned long long num_combs = combinations(n, i);
#pragma omp parallel for shared(exit, flag, out_x)
        for (unsigned long long cid = 0; cid < num_combs; cid++) {
            if (flag || exit)
                continue;
            vector<int> v_indices = ith_combination(n, i, cid);

            // Check if n is even and if we are at the half-way point.
            if (n % 2 == 0 && i == w && v_indices[0] != 0) {
                flag = true;
                continue;
            }

            vector<uint8_t> x(n, 0);
            for (int idx : v_indices) {
                if (idx < 0 || idx >= n) {
                    cerr << "Error: índice fuera de rango " << idx
                         << " para n=" << n << ", i=" << i << ", cid=" << cid
                         << endl;
                    abort();
                }
                x[idx] = 1;
            }

            if (i < min_len) {
                vector<uint8_t> one_minus_x(n);
                for (int k = 0; k < n; ++k) {
                    one_minus_x[k] = 1 - x[k];
                }
                int w_val = n - i;
                if (evaluate_f(one_minus_x, w_val) == 0) {
#pragma omp critical(outXCopy)
                    if (!exit) {
                        exit = true;
                        copy(x.begin(), x.end(), out_x.begin());
                    }
                }
            } else {
                if (evaluate_f(x, i) == 0) {
                    // if (evaluate_f_old(x) == 0) {
                    vector<uint8_t> one_minus_x(n);
                    for (int k = 0; k < n; ++k) {
                        one_minus_x[k] = 1 - x[k];
                    }
                    int w_val = n - i;
                    if (evaluate_f(one_minus_x, w_val) == 0) {
                        // if (evaluate_f_old(one_minus_x) == 0) {
#pragma omp critical(outXCopy)
                        if (!exit) {
                            exit = true;
                            copy(x.begin(), x.end(), out_x.begin());
                        }
                    }
                }
            }
        }
        if (exit)
            return {false, out_x};
    }
    return {true, {}}; // Return true and an empty vector if no such x is found
}

// Returns a vector of vectors (hyperedges) sorted by their elements
vector<vector<int>> hypergraph::get_ordered_list_psi() {
    vector<vector<int>> ordered_psi(psi.begin(), psi.end());
    // The custom comparator `VectorIntCompare` already sorts by size then
    // lexicographically So simply copying to vector is sufficient.
    return ordered_psi;
}

// Returns the hyperedge with the minimum length
vector<int> hypergraph::min_length_hyperedge() {
    if (psi.empty()) {
        return {};
    }
    return *min_element(psi.begin(), psi.end(),
                        [](const vector<int> &a, const vector<int> &b) {
                            return a.size() < b.size();
                        });
}

// Returns the hyperedge with the maximum length
vector<int> hypergraph::max_length_hyperedge() {
    if (psi.empty()) {
        return {};
    }
    return *max_element(psi.begin(), psi.end(),
                        [](const vector<int> &a, const vector<int> &b) {
                            return a.size() < b.size();
                        });
}

void hypergraph::save(const string &fname) {
    json d;
    d["n"] = n;
    d["psi"] = psi;

    ofstream fileOutput(fname);
    if (fileOutput.is_open()) {
        fileOutput << d.dump();
        fileOutput.close();
        cout << "Data succesfully saved in " << fname << endl;
    } else {
        cerr << "Error: Impossibile to open the file " << fname
             << " for writing." << endl;
    }
}

void hypergraph::load(const string &fname) {
    /*
        Load the hypergraph from a JSON file.
        The JSON file should contain the number of vertices 'n' and the set of
       hyperedges 'psi'. The 'psi' should be a list of lists, where each inner
       list represents a hyperedge.
    */
    if (fname.size() >= 4 && fname.substr(fname.size() - 4) == ".dat") {
        load_uno_dat(fname);
        return;
    }
    ifstream fileInput(fname);
    if (fileInput.is_open()) {
        json d;
        try {
            fileInput >> d; // Parsifica il JSON dal file
            fileInput.close();
            n = d["n"].get<int>();
            psi.clear();
            for (const auto &e : d["psi"]) {
                vector<int> edge = e.get<vector<int>>();
                psi.insert(edge);
            }
        } catch (const json::parse_error &e) {
            cerr << "Error parsing JSON: " << e.what() << endl;
            return;
        } catch (const json::type_error &e) {
            cerr << "Errore of JSON type: " << e.what() << endl;
            return;
        }
    } else {
        cerr << "Error: Cannot open file " << fname << " for reading." << endl;
        return;
    }
}

void hypergraph::load_uno_dat(const string &filename) {

    ifstream infile(filename);
    if (!infile) {
        cerr << "Error: Cannot open file " << filename << endl;
        return;
    }
    set<int> vertices;
    string line;
    psi.clear();
    while (getline(infile, line)) {
        istringstream iss(line);
        vector<int> edge;
        int val;
        while (iss >> val) {
            edge.push_back(val);
            vertices.insert(val); // Collect unique vertices
        }
        if (!edge.empty()) {
            psi.insert(edge);
        }
    }
    n = vertices.size(); // Set n to the number of unique vertices
    infile.close();
    cout << "Loaded " << psi.size() << " hyperedges with " << n
         << " unique vertices from " << filename << endl;
}

map<string, double> hypergraph::stat() {
    map<string, double> st;
    st["n"] = static_cast<double>(n);
    st["m"] = static_cast<double>(psi.size());

    if (!psi.empty()) {
        st["min_edge"] = static_cast<double>(min_length_hyperedge().size());
        st["max_edge"] = static_cast<double>(max_length_hyperedge().size());

        double sum_len = 0;
        for (const auto &e : psi) {
            sum_len += e.size();
        }
        st["avg"] = sum_len / st["m"];
    } else {
        st["min_edge"] = 0.0;
        st["max_edge"] = 0.0;
        st["avg"] = 0.0;
    }
    return st;
}

// ---
// Global function reduce_minimal
// ---

void reduce_minimal(vector<vector<int>> &E) {
    set<vector<int>, VectorIntCompare> edges_to_remove;

    for (int i = 0; i < E.size(); ++i) {
        set<int> a = vec_to_set(E[i]);
        for (int j = i + 1; j < E.size(); ++j) {
            set<int> b = vec_to_set(E[j]);

            bool a_is_subset_b = true;
            for (int val : a) {
                if (b.find(val) == b.end()) {
                    a_is_subset_b = false;
                    break;
                }
            }

            bool b_is_subset_a = true;
            for (int val : b) {
                if (a.find(val) == a.end()) {
                    b_is_subset_a = false;
                    break;
                }
            }

            if (a_is_subset_b) {
                edges_to_remove.insert(E[j]);
            } else if (b_is_subset_a && a != b) {
                edges_to_remove.insert(E[i]);
            }
        }
    }

    // Remove elements from E
    E.erase(remove_if(E.begin(), E.end(),
                      [&](const vector<int> &edge) {
                          return edges_to_remove.count(edge);
                      }),
            E.end());
}

// ---
// Global function algorithm_A
// ---

bool algorithm_A(vector<vector<int>> F, vector<vector<int>> G) {
    /*
        Implementation of the algorithm A of paper
        [1] Fredman, M.L., Khachiyan, L.: On the complexity of dualization of
        monotone disjunctive normal forms. J. Algorithms 21(3), 618–628 (1996)

        Returns
        -------
        bool.
            true if G is the dual of F
    */
    // Deal with the base case
    if (F.size() < 2 && G.size() < 2) {
        if (F.size() == 0 && G.size() == 0) {
            return true;
        } else if (F[0] == G[0]) {
            return true;
        } else {
            return false;
        }
    }

    reduce_minimal(F);
    reduce_minimal(G);

    // Check condition (1.2) of the paper
    map<int, int> V_F;
    map<int, int> V_G;
    int x_i = 0, m_i = 0;
    int x_j = 0, m_j = 0;
    int I_m = 0, J_m = 0;

    for (const auto &e : F) {
        if (e.size() > I_m) {
            I_m = e.size();
        }
        for (int v : e) {
            V_F[v]++;
            if (V_F[v] > m_i) {
                x_i = v;
                m_i = V_F[v];
            }
        }
    }
    for (const auto &e : G) {
        if (e.size() > J_m) {
            J_m = e.size();
        }
        for (int v : e) {
            V_G[v]++;
            if (V_G[v] > m_j) {
                x_j = v;
                m_j = V_G[v];
            }
        }
    }

    // Check if key sets are equal
    if (V_F.size() !=
        V_G.size()) { // If their sizes are different, they cannot be equal
        return false;
    }
    for (const auto &pair :
         V_F) { // They have equal size, so check if all keys in V_F are in V_G
        if (V_G.find(pair.first) == V_G.end()) {
            return false;
        }
    }

    // Check condition (1.3)
    if (I_m > G.size() || J_m > F.size()) {
        return false;
    }

    /**********************************
     *  Check condition (2.1)
     ***********************************
     */
    /* old code using double data type
    double E_sum = 0.0;
    for (const auto& e : F) {
        E_sum += pow(2, -static_cast<double>(e.size()));
    }
    for (const auto& e : G) {
        E_sum += pow(2, -static_cast<double>(e.size()));
    }
    if (E_sum < 1.0) {
        return false;
    }
    */

    mpz_class E_sum = 0;
    mpz_class M, E;
    int n = V_F.size(); // Assuming n is the number of unique vertices
    mpz_ui_pow_ui(M.get_mpz_t(), 2, n); // M = 2^n

    for (const auto &e : F) {
        mpz_ui_pow_ui(E.get_mpz_t(), 2, n - e.size()); // E = 2^(n - |e|)
        E_sum += E;
    }
    for (const auto &e : G) {
        mpz_ui_pow_ui(E.get_mpz_t(), 2, n - e.size()); // E = 2^(n - |e|)
        E_sum += E;
    }

    if (E_sum < M) { // If the sum of the weights is less than 2^n, return false
        return false;
    }

    // Detrmines the vertex of maximum frequency in F and G
    if (m_i < m_j) {
        x_i = x_j;
    }

    // Create subproblems F_0, F_1, G_0, G_1
    vector<vector<int>> F_0;
    vector<vector<int>> F_1;
    for (const auto &e : F) {
        vector<int> e_copy = e;
        auto it = remove(e_copy.begin(), e_copy.end(), x_i);
        bool found_x_i = (it != e_copy.end());

        if (found_x_i) {
            e_copy.erase(it, e_copy.end());
            if (!e_copy.empty()) {
                F_0.push_back(e_copy);
            }
        } else {
            F_1.push_back(e_copy);
        }
    }

    vector<vector<int>> G_0;
    vector<vector<int>> G_1;
    for (const auto &e : G) {
        vector<int> e_copy = e;
        auto it = remove(e_copy.begin(), e_copy.end(), x_i);
        bool found_x_i = (it != e_copy.end());

        if (found_x_i) {
            e_copy.erase(it, e_copy.end());
            if (!e_copy.empty()) {
                G_0.push_back(e_copy);
            }
        } else {
            G_1.push_back(e_copy);
        }
    }

    vector<vector<int>> G_0G_1;
    if (!G_0.empty()) {
        G_0G_1.insert(G_0G_1.end(), G_0.begin(), G_0.end());
        if (!G_1.empty()) {
            G_0G_1.insert(G_0G_1.end(), G_1.begin(), G_1.end());
        }
    }

    vector<vector<int>> F_0F_1;
    if (!F_0.empty()) {
        F_0F_1.insert(F_0F_1.end(), F_0.begin(), F_0.end());
        if (!F_1.empty()) {
            F_0F_1.insert(F_0F_1.end(), F_1.begin(), F_1.end());
        }
    }

    // Recursive calls
    if (!algorithm_A(F_1, G_0G_1)) {
        return false;
    }
    return algorithm_A(G_1, F_0F_1);
}

void print_vector_of_vectors(const vector<vector<int>> &F, string l) {
    cout << "Vector: " << l << endl;
    for (const auto &edge : F) {
        cout << "{ ";
        for (int v : edge) {
            cout << v << " ";
        }
        cout << "}" << endl;
    }
}

// utility function for date and time string
string get_datetime_string() {
    auto now = chrono::system_clock::now();
    time_t now_time = chrono::system_clock::to_time_t(now);
    tm *now_tm = localtime(&now_time);
    ostringstream oss;
    oss << put_time(now_tm, "%Y_%m_%d_%H_%M_%S");
    return oss.str();
}
