# dualization
Implementation in c++ of algorithms for solving the dualization problem for positive Boolean functions

# Description
The code provided is the implementation af some algorithms for the decision version of the dualization problem for positive Boolean function. There are tree classes: "bipartite_graph" used to deal with the bipartite representation of a hypergraph. In this class are implemented the algorithm:
* "compute_number_hitting_sets" which computes the number of hitting sets of a hypergraph
  
The class "hypergraph" is used to store the raw hypergraph as a set of integer vectors. In this class are implemented the algorithms:
* "generate_random_comb_new" generate hypergraph with specified type of hyperedges using uniform random distribution
* "generate_random_exponential" generate hypergraph with specified type of hyperedges using exponential random distribution
* "search_x" search for a $x$ such that $f(x)=f(\overline{x})$
* "sum_f" computes $\sum_{x=0}^{2^{n}-1} [1-f(x)]$ 

The class "bipartite_02" which is derived from class bipartite_graph and implements the algorithm
* "find_uncovering_hitting_set" which find an uncovering hitting set if any exists in the hypergraph
  
In the file "bipartite_hypergraph.cpp" there is also an implementation of the algorithm A of paper [1]    

The directory "saved_hypergraphs" contains several hypergraphs used in the experiments. The JSON file is a simple maps containing two keys: "n": the number of vertices of the hypergraph; "psi": a vector of integer vectors.

# Dependencies
The GMP libreries are used for arbitrary aritmetic precision. The header "json.hpp2 is used for dealing with JSON file format.

# Usage
Run "experiment_prod_ver1.5" to execute some algorithms on some hypergraphs.
Run "generate_hypergraph_ver1.0" to randomly generate and save hypergraphs.

A Makefile is provided to build the code.


# Bibliography
[1] Fredman, M.L., Khachiyan, L.: On the complexity of dualization of 
      monotone disjunctive normal forms. J. Algorithms 21(3), 618–628 (1996)
