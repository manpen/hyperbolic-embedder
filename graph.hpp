#pragma once

#include <string>
#include <unordered_map>
#include <vector>

using std::vector;
using std::ostream;
using std::string;
using std::unordered_map;

class Graph
{
 public:
  int n;
  vector<vector<int>> edges;
  vector<std::string> labels;
  
  Graph(int n) {
    this->n = n;
    //edges.resize(n);
    //for(int i = 0; i < n; ++i) labels.push_back(std::to_string(i));
  }
  
  Graph(const Graph& G);
  virtual ~Graph(){}
  static Graph fromFile(const string& filename,
                        unordered_map<std::string, int>* label_to_node);

  // Changes the number of nodes. Adds zero-degree vertices, if new n is larger
  // than old n. Deletes the last vertices, if new n is smaller then old n.
  // Warning/TODO: Does not delete edges to deleted nodes.
  void resize(int n);

  Graph subgraph(vector<bool> vertices);
  // Returns the subgraph of the largest connected component in the graph.
  Graph giantSubgraph();
  // Returns a simple subgraph, i.e. no self-loops and no double edges.
  Graph simpleSubgraph();
  
  double averageDegree() const;
  vector<int> degHisto() const;

  vector<int> components() const;
  vector<bool> giant(vector<int> components) const;
  vector<int> smallComponentHisto(vector<int> components, vector<bool> giant);
  int sizeOfGiant(vector<bool> giant) const;
  
  bool adjacent(int u, int v) const;

  // Permutes nodes according to the given permutation. If permutation[i] == -1,
  // discards node i.
  void permuteGraph(const vector<int>& permutation);
  // Subsamples the graph's edges, keeping each edge with probability p.
  Graph subsample(double p);

  // Reorders the adjacency lists s.t. they are ordered according to the given
  // permutation of nodes. Does not change the graph structure.
  // E.g. if perm is [4,3,5,1,2] then a node with neighbors 1,2,3 stores
  // the neighbors as [3,1,2].
  void ReorderEdges(const vector<int>& permutation);

  // Reorders the nodes by decreasing degree. edges[0] is the node with the
  // largest degree.
  // new_permutation will contain the reordering of the nodes. permutation[i]=j
  // means that old node i is now node j.
  virtual void sortByDegrees(vector<int>* new_permutation);

  // Saves this graph to a file with the given name. Saves each edge only once,
  // i.e. considers the graph to be undirected.
  // Ignores all edges (u,v) where either ignore_nodes[u] == true or
  // ignore_nodes[v] == true.
  virtual void printToFile(const char* filename,
                           const vector<bool>* ignore_nodes = NULL,
                           bool use_labels = true);

  // Counts the number of common neighbors of u and v
  int commonNeighbors(int u, int v);

 private:
  //void mark(int i, int T, vector<int> & comp);
};
