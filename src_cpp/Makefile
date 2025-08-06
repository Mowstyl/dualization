CXX = g++
CXXFLAGS = -O3 -std=c++17 -march=native -flto -Wall -Wextra
LDFLAGS = -flto -lgmp -lgmpxx

TARGETS = generate_hypergraph_ver1.0 experiment_prod_ver1.5

all: $(TARGETS)

generate_hypergraph_ver1.0: generate_hypergraph_ver1.0.o bipartite_hypergraph.o bipartite_02.o
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDFLAGS)

experiment_prod_ver1.5: experiment_prod_ver1.5.o bipartite_hypergraph.o bipartite_02.o
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDFLAGS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f *.o $(TARGETS)

run_generate:
	./generate_hypergraph_ver1.0

run_experiment:
	./experiment_prod_ver1.5
