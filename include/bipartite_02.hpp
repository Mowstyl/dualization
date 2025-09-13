#include "bipartite_hypergraph.hpp"

class bipartite_02 : public bipartite_graph {
  public:
    vector<int> counter_example; // counter example
    int min_level;
    mpz_class count;
    const vector<vector<int>> *Edges_copy;
    bipartite_02(const vector<vector<int>> &H);
    bipartite_02(const bipartite_02 &other);
    bipartite_02 &operator=(const bipartite_02 &other);
    bool is_uncovering(const vector<int> &ce) const;
    tuple<map<int, vector<int>>, map<int, vector<int>>, set<int>>
    compute_H_2(const vector<int> &p, int e_id);
    bool find_uncovering_hitting_set(int level);
};
