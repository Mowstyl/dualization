#include "bipartite_02.hpp"

bipartite_02::bipartite_02(const vector<vector<int>> &H)
    : bipartite_graph(H), Edges_copy(&H), counter_example() {
    // Constructor that initializes the bipartite graph with the given
    // hypergraph Edges_copy is initialized in the initializer list
}
// This constructor initializes the bipartite graph with the hypergraph edges
// and can be used to create a bipartite_02 object from a hypergraph
// representation This is useful for algorithms that require a bipartite
// representation of the hypergraph

bipartite_02::bipartite_02(const bipartite_02 &other)
    : bipartite_graph(other), Edges_copy(other.Edges_copy) {
    // Copy constructor that initializes the bipartite graph from another
    // bipartite_02 object
    counter_example = other.counter_example;
}

bipartite_02 &bipartite_02::operator=(const bipartite_02 &other) {
    if (this != &other) {
        bipartite_graph::operator=(other);
        counter_example = other.counter_example;
        Edges_copy = other.Edges_copy;
    }
    return *this;
}

bool bipartite_02::is_uncovering(const vector<int> &ce) const {
    /*
    Checks if the counter example is uncovering. That is if there is no edge
    that is contained in the counter example. To work correctly, the counter
    example must be sorted. Also the set of hyperedges must be sorted by length
    and hyperedges of the same length must be sorted by their vertices.

        Parameters
    ----------
    ce : list or tuple
        the set of vertices to be checked.
    Returns
    -------
    bool
        True if the set is uncovering, False otherwise.
    */
    size_t ls = ce.size();
    // sort the counter example
    vector<int> ce_sorted = ce;
    sort(ce_sorted.begin(), ce_sorted.end());

    for (vector<int> e : *Edges_copy) {
        /* Since the hyperedges in Edges_copy are sorted by length, if the
          counter example is smaller than the hyperedge e, e and all othere
          subsequent hyperedges cannot be contained in it
        */
        if (ls < e.size()) {
            return true;
        }

        size_t i = 0, j = 0, count = 0;
        while (i < ls && j < e.size()) {
            if (ce_sorted[i] < e[j]) {
                ++i;
            } else if (ce_sorted[i] > e[j]) {
                ++j;
            } else {
                ++count;
                ++i;
                ++j;
            }
        }
        if (count == e.size()) {
            return false;
        }
    }
    return true;
}

static vector<int> set_difference(const vector<int> &a, const vector<int> &b) {
    vector<int> result;
    set<int> b_set(b.begin(), b.end());
    for (int x : a) {
        if (b_set.find(x) == b_set.end()) {
            result.push_back(x);
        }
    }
    return result;
}

tuple<map<int, vector<int>>, map<int, vector<int>>, set<int>>
bipartite_02::compute_H_2(const vector<int> &p, int e_id) {
    /*
            Parameters
        ----------
        p : vector of integer
            the set of vertices neede to be removed .
        e_id : int
            index of the edges considered.
        Computes $H_2$ as reported in the paper. First $H_1=H \setmins p$
        is computed. Then $H_2 = H_1 \setminus N_{H_1}(s)$ is computed, where
       $s$ is the edge with index e_id. The function returns the backup of the
       edges, adjacency lists and vertices before the deletion.

        Returns
        -------

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

    /*  Delete the vertices in V(H_1) \setminus (V(H_2) \cup E_copy) */
    vector<int> deleted_vertices_isolated;

    for (int v : Vertices) {
        if (!Adj.count(v) ||
            Adj[v].empty()) { // Check if adjacency list exists and is empty
            deleted_vertices_isolated.push_back(v);
        }
    }

    for (int v : deleted_vertices_isolated) {
        Vertices.erase(v);
        Adj.erase(v);
    }
    /*
    print_graph("Bipartite Graph after compute_H_2: \n ");
    cout << "Edges size after compute H_2: " << Edges.size() << endl;
    */
    return make_tuple(edges_backup, adj_backup, vertices_backup);
}

bool bipartite_02::find_uncovering_hitting_set(int level) {

    if (Edges.empty()) {
        // print_vector(counter_example, "Counter example 0: ");
        return is_uncovering(counter_example);
    }
    if (Edges.size() == 1) {
        // if there is only one edge, we can use only one vertex at a time from
        // it
        for (int v : Edges.begin()->second) {
            counter_example.push_back(v);
            // print_vector(counter_example, "Counter example 1: ");

            if (is_uncovering(counter_example)) {

                return true; // Found a minimal hitting set
            }

            counter_example.pop_back(); // Backtrack
        }
        return false;
    }

    int min_len = numeric_limits<int>::max();
    int e_min_len_id = -1; // The id of the edge with minimum cardinality

    for (const auto &pair_e : Edges) {
        if (pair_e.second.size() < min_len) {
            e_min_len_id = pair_e.first;
            min_len = pair_e.second.size();
        }
    }

    vector<int> E_current_edge = Edges[e_min_len_id];
    vector<int> elements_for_combinations(E_current_edge.begin(),
                                          E_current_edge.end());
    vector<int> p;
    vector<int> e_p;

    int min_iter = min_len;
    // min_iter = 1;
    /*
    if (level == 0) {
        min_iter = min_len-1;
    } else {
        min_iter = min_len - min_level; // Use the minimum level found in the
    previous recursion if (min_iter <0) { min_iter = 0; // Ensure we don't have
    negative iterations
        }
    }
    */

    for (int i = 1; i <= min_iter; ++i) {
        for (CombinationGenerator<int> gen(elements_for_combinations, i);
             !gen.is_done();) {
            e_p = gen.next();
            /*
            if (level == 0) {
                min_level = j; // Set the minimum level for the successive
            recursion levels
            }
            */

            p = set_difference(E_current_edge, e_p);
            counter_example.insert(counter_example.end(), e_p.begin(),
                                   e_p.end());

            // print_vector(p, "p: ");
            // print_vector(counter_example, "counter_example: ");

            auto [e_b, a_b, v_b] = this->compute_H_2(p, e_min_len_id);

            bool n_r = this->find_uncovering_hitting_set(level + 1);

            if (n_r)
                return true; // Found an uncovering hitting set

            // Revert changes back to *this for the next iteration
            // Restore the edges, adjacency lists, and vertices

            for (const auto &pair : e_b) {
                this->Edges[pair.first] = pair.second;
            }
            for (const auto &pair : a_b) {
                this->Adj[pair.first] = pair.second;
            }
            for (int v : v_b) {
                this->Vertices.insert(v);
            }
            // Restore the state of the counter_example
            // by removing the last elements that were added
            counter_example.erase(counter_example.end() - e_p.size(),
                                  counter_example.end());
        }
    }
    return false;
}

/*
int main() {
    test_bipartite_02();  return 0;
}

*/