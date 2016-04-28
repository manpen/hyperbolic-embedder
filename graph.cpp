#include "graph.hpp"
#include "random.hpp"

#include <algorithm>
#include <iostream>
#include <fstream>
#include <map>
#include <queue>
#include <stack>
#include <glog/logging.h>
#include <cstdio>
#include <cmath>
#include <set>
#include <unordered_set>

#define FOR(i,n) for(int (i)=0;(i)<(n);(i)++)
#define FORB(i,a,n) for(int (i)=(a);(i)<(n);(i)++)
#define FOREACH(it,c) for(__typeof((c).begin()) it=(c).begin();it!=(c).end();++it)
#define PB push_back
#define MP make_pair

using std::min;
using std::max;
using std::set;
using std::map;
using std::queue;
using std::stack;
using std::swap;
using std::endl;
using std::cout;

Graph::Graph(const Graph& G) {
  n = G.n;
  edges.resize(n);
  FOR(i,n)
    FOR(j,G.edges[i].size())
      edges[i].PB(G.edges[i][j]);
}

Graph Graph::fromFile(const string& filename,
                      unordered_map<string, int>* label_to_node) {
  Graph g(0);
  std::ifstream file;
  file.open(filename.c_str());
  CHECK(file.good()) << "Maybe the file " << filename << " doesn't exist?";
  label_to_node->clear();
  g.labels.clear();

  // read all nodes
  int current_node = 0;
  string current_label;
  while (!file.eof()) {
    if ((file >> std::ws).peek() == '#') { // ignore rest of the line
      std::getline(file, current_label);
      continue;
    }
    file >> current_label;
    // ignore directory structure
    current_label = current_label.substr(0, current_label.find("/"));
    if (label_to_node->count(current_label) == 0) {
      // new node
      (*label_to_node)[current_label] = current_node;
      ++current_node;
      g.labels.PB(current_label);
    }
  }
  file.close();

  // now construct graph
  g.resize(current_node);

  file.open(filename.c_str());
  string u, v;
  while (!file.eof()) {
    if ((file >> std::ws).peek() == '#') { // ignore rest of the line
      std::getline(file, current_label);
      continue;
    }
    file >> u; file >> v;
    u = u.substr(0, u.find("/"));
    v = v.substr(0, v.find("/"));
    if ((*label_to_node).find(u) != label_to_node->end() &&
        (*label_to_node).find(v) != label_to_node->end()) {
      g.edges[(*label_to_node)[u]].PB((*label_to_node)[v]);
      g.edges[(*label_to_node)[v]].PB((*label_to_node)[u]);
    }
  }
  file.close();

  return g;
}


void Graph::resize(int n) {
  int old_n = this->n;
  this->n = n;
  edges.resize(n);
  labels.resize(n,"");

  FORB(i,old_n,n) if (labels[i].empty()) labels[i] = std::to_string(i);
}

Graph Graph::subgraph(vector<bool> vertices) {
  vector<int> index (n+1);
  index[0] = 0;
  FORB(i,1,n+1) index[i] = index[i-1] + (vertices[i-1] ? 1 : 0);
  int m = index[n];
  Graph G (m);
  FOR(i,n){
    if (vertices[i]) {
      G.labels[index[i]] = labels[i];
      FOR(j,(int)edges[i].size()) {
        int jj = edges[i][j];
        if (vertices[jj])
          G.edges[index[i]].PB(index[jj]);
      }
    }
  }
  return G;
}


double Graph::averageDegree() const {
  int m = 0;
  FOR(i,n) m += (int) edges[i].size();
  return ((double)m) / ((double)n);
}

vector<int> Graph::degHisto() const {
  int maxDeg = 0;
  FOR(i,n) maxDeg = max(maxDeg, (int) edges[i].size());
  vector<int> histo (maxDeg+1);
  FOR(i,n) histo[edges[i].size()]++;
  return histo;
}

vector<int> Graph::components() const {
  vector<int> comp(n,-1);
  int T = 0;
  FOR(i,n) 
    if (comp[i] == -1) {
      stack<int> S;
      S.push(i);
      while(!S.empty()) {
        int v = S.top(); S.pop();
        if (comp[v] == -1) {
          comp[v] = T;
          FOREACH(it,edges[v])
            S.push(*it);
        }
      }
      
      T++;
    }
  return comp;
}

vector<bool> Graph::giant(vector<int> components) const {
  vector<int> size (n,0);
  FOR(i,n) size[components[i]]++;
  int maxcomp = 0;
  FOR(i,n) 
    if (size[i] > size[maxcomp])
      maxcomp = i;
  vector<bool> giant (n);
  FOR(i,n)
    giant[i] = (components[i] == maxcomp);
  return giant;
}

vector<int> Graph::smallComponentHisto(vector<int> components, vector<bool> giant) {
  vector<int> size (n,0);
  int maxsize = 0;
  FOR(i,n) 
    if (!giant[i]) {
      size[components[i]]++;
      maxsize = max(maxsize, size[components[i]]);
    }
  vector<int> histo (maxsize+1,0);
  FOR(i,n) histo[size[i]]++;
  histo[0] = 0;
  return histo;
}

int Graph::sizeOfGiant(vector<bool> giant) const {
  int size = 0;
  FOR(i,n) size += (giant[i]) ? 1 : 0;
  return size;
}

bool Graph::adjacent(int u, int v) const {
  if (edges[u].size() > edges[v].size()) return adjacent(v,u);
  FOREACH(w,edges[u])
    if (*w == v) return true;
  return false;
}
Graph Graph::giantSubgraph() {
  vector<int> comp = components();
  vector<bool> gcc = giant(comp);
  return subgraph(gcc);
}

Graph Graph::simpleSubgraph() {
  vector<bool> seen (n,false);
  Graph G (n);
  FOR(i,n) {
    G.labels[i] = labels[i];
    for (int node : edges[i])
      if (!seen[node] && (node != i)) {
        seen[node] = true;
        G.edges[i].PB(node);
      }
    std::fill(seen.begin(), seen.end(), false);
  }
  return G;
}

void Graph::permuteGraph(const vector<int>& permutation) {
  int num_subsampled_nodes = 0;
  FOR(i,permutation.size())
    if (permutation[i] >= 0)
      ++num_subsampled_nodes;

  vector<vector<int>> new_edges(num_subsampled_nodes);
  vector<std::string> new_labels(num_subsampled_nodes);
  FOR(i,n) {
    if (permutation[i] >= 0) {
      new_edges[permutation[i]].swap(edges[i]);
      new_labels[permutation[i]] = labels[i];
    }
  }


  FOR(i,num_subsampled_nodes) {
    vector<int> new_edges_i;
    FOR(j,new_edges[i].size())
      if (permutation[new_edges[i][j]] >= 0)
        new_edges_i.PB(permutation[new_edges[i][j]]);
    new_edges[i].swap(new_edges_i);
  }

  edges.swap(new_edges);
  labels.swap(new_labels);
  n = (int)edges.size();
}

Graph Graph::subsample(double p) {
  Graph H(n);

  FOR(i,n) {
    FOR(j,edges[i].size()) {
      if (edges[i][j] <= i)
        continue; // Saw this edge already
      if (randdblpos() < p) {
        H.edges[i].PB(edges[i][j]);
        H.edges[edges[i][j]].PB(i);
      }
    }
  }

  return H;
}

void Graph::sortByDegrees(vector<int>* new_permutation) {
  // calculate the histogram of key frequencies:
  vector<int> count(n,0);
  FOR(i,n)
    count[edges[i].size()]++;

  // calculate the starting index for each key:
  int total = 0;
  for (int i = n-1; i >= 0; --i) {
    int oldCount = count[i];
    count[i] = total;
    total += oldCount;
  }

  new_permutation->clear();
  new_permutation->resize(n);

  // copy to output array, preserving order of inputs with equal keys:
  FOR(i,n) {
    int deg = (int)edges[i].size();
    (*new_permutation)[i] = count[deg];
    count[deg] += 1;
  }

  permuteGraph(*new_permutation);
}

void Graph::printToFile(const char* filename,
                        const vector<bool>* ignore_nodes,
                        bool use_labels) {
    std::ofstream file;
    file.open(filename);

    FOR(i,n) {
      for (int j : edges[i]) {
        if (ignore_nodes == NULL ||
            (!ignore_nodes->at(i) && !ignore_nodes->at(j)))
          if (j >= i) {
            if (use_labels) file << labels[i] << "\t" << labels[j] << std::endl;
            else            file <<        i  << "\t" <<        j  << std::endl;
          }
      }
    }

    file.close();
}

int Graph::commonNeighbors(int u, int v) {
  int common_neighbors = 0;
  FOR(i,edges[u].size())
    FOR(j,edges[v].size())
      if (edges[u][i] == edges[v][j])
        ++common_neighbors;

  return common_neighbors;
}

void Graph::ReorderEdges(const vector<int>& permutation) {
  vector<vector<int>> new_edges(n);

  FOR(i,n) {
    FOR(j,edges[permutation[i]].size()) {
      new_edges[edges[permutation[i]][j]].PB(permutation[i]);
    }
  }

  edges.swap(new_edges);
}

namespace{
struct lexcomp {
  bool operator() (const vector<int>& lhs, const vector<int>& rhs) const {
    int i = 0;
    while (i < rhs.size()) {
      if (i >= lhs.size() || rhs[i]<lhs[i]) return false;
      else if (lhs[i]<rhs[i]) return true;
      ++i;
    }
    return (i != lhs.size());
  }
};
}