// singleexperiment.cpp
// 
// Used for running some experiments on random graphs. Usage is:
// ./singleexperiment --rgmodel=chunglu --average_degree
#include <algorithm>
#include <cmath>
#include <ctime>
#include <gflags/gflags.h>
#include <glog/logging.h>
#include <iostream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "graph.hpp"
#include "hyperbolic.hpp"
#include "hyperboliclinear.hpp"
#include "random.hpp"
#include "springembedder.hpp"
#include "NLEHelper.hpp"

#define FOR(i,n) for(int (i)=0;(i)<(n);(i)++)
#define FORB(i,a,n) for(int (i)=(a);(i)<(n);(i)++)
#define FOREACH(it,c) for(__typeof((c).begin()) it=(c).begin();it!=(c).end();++it)
#define PB push_back
#define MP make_pair

using std::unordered_map;
using std::unordered_set;
using std::cout;
using std::endl;
using std::vector;
using std::string;

DEFINE_int64(seed, 0, "Seed for RNG. If 0, uses the current time of day.");
DEFINE_string(input, "", "Graph to be embedded. Expects a file with an edge "
              "list as input. "
              "I.e. on each line it reads the two endpoints of an edge. By "
              "design, does not support nodes with no neighbors. Note that "
              "only the giant component of the given graph is embedded.");

DEFINE_string(generate, "", "If a non-empty string is given, will generate a "
              "hyperbolic random graph for embedding instead of reading a "
              "graph from an input file. For reference, the generated graph "
              "is stored in name-links.txt, and the locations are stored in "
              "name-coords.txt, where name is the supplied string.");

DEFINE_int32(n, 1000, "Number of vertices for hyperbolic random graph");
DEFINE_double(C, -1, "Adjusts the average degree of the generated hyperbolic "
              "random graph, can be any real number. Radius of the disc will "
              "be R=2 * log(n) + C. Large C means small average degree.");
DEFINE_double(T, 0.1,
    "Adjusts probability that edge exists depending on distance of two nodes.");
DEFINE_double(alpha, 0.75, "Adjusts power-law exponent of the generated "
              "hyperbolic random graph.");

DEFINE_string(embed, "embedded", "Embeds the given graph and saves the "
              "embedded graph and the inferred coordinates to name-links.txt "
              "and name-coordinates.txt, where name is the supplied string.");
DEFINE_string(springembed, "", "If not empty, (also) embeds the graph using "
              "the spring embedder and saves the coordinates to the given "
              "location.");


int main(int argc, char* argv[]) {
  google::ParseCommandLineFlags(&argc, &argv, true);
  google::InitGoogleLogging(argv[0]);
  if (FLAGS_seed == 0)
    FLAGS_seed = std::time(NULL);
  initRNG(FLAGS_seed);
  srand((unsigned int) FLAGS_seed);

  cout << "Seed: " << FLAGS_seed << endl;

  Graph G(0);

  if (!FLAGS_generate.empty()) {
    double R = 2 * log(FLAGS_n) + FLAGS_C;
    HyperbolicLinear H = HyperbolicLinear::linearSampling(FLAGS_n, R, FLAGS_alpha, FLAGS_T);
    H.printToFile(FLAGS_generate.c_str());
    G = H.simpleSubgraph();
  } else if (!FLAGS_input.empty()) {
    unordered_map<std::string, int> label_to_node;
    G = Graph::fromFile(FLAGS_input, &label_to_node);
    G = G.simpleSubgraph();
  } else {
    cout << "Please either specify an input graph using --input or use "
         << "--generate to generate a random hyperbolic graph. Call with "
         << "--help for more information." << endl;
    exit(0);
  }

  // Compute giant and sort by degrees
  G = G.giantSubgraph();
  vector<int> perm;
  G.sortByDegrees(&perm);

  cout << "Number of nodes in giant: " << G.n << endl;
  cout << "Avg degree: " << G.averageDegree() << endl;

  double R, T, alpha;
  double n_orig, m_orig;
  vector<double> est_r(G.n);
  NLEHelper::estimateHyperbolicParameters(G, &T, &n_orig, &m_orig, &alpha, &R,
                                          &est_r);

  HyperbolicLinear HE(G.n,R,alpha,T);
  HE.edges = G.edges;
  HE.labels = G.labels;
  HE.pts.resize(G.n);
  FOR(i,HE.n) {
    HE.pts[i].r = est_r[i];
    HE.pts[i].phi = M_PI * 1.5;
  }

  if (!FLAGS_embed.empty() && FLAGS_springembed.empty()) {
    cout << "Embedding graph in the hyperbolic plane..." << endl;
    HE.hyperbolicEmbedding();
    HE.printToFile(FLAGS_embed.c_str());
    cout << "Embedded Log-likelihood: " << HE.computeEnergy(0) << endl;
  }
  if (!FLAGS_springembed.empty()) {
    cout << "Embedding graph in the hyperbolic plane using a spring embedder..."
         << endl;
    SpringEmbedder SE(HE);
    SE.run(1000, 1, 0, true);
    HE.printToFile(FLAGS_springembed.c_str());
  }

  cout << "End." << endl;
}

