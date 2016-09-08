#pragma once

#include "graph.hpp"
#include <cmath>
#include <set>
#include <string>
#include <ostream>


class HYPT {
 public:
  double r, phi;
  
  HYPT(double r, double phi);
    
  HYPT() : r(0), phi(0) {};

  static HYPT randomHYPT(double R, double alpha);
  static double dist(HYPT p, HYPT q);
  static double dist(double r1, double r2, double phi1, double phi2);
};

ostream& operator<<(ostream& os, const HYPT& pt);

class Hyperbolic : public Graph {
 public:
  vector<HYPT> pts;
  double R;
  double alpha;
  double T;

  static double prob(double dist, double R, double T) {
    if (T == 0.) return (dist <= R) ? 1. : 0.;
    return 1./ ( 1. + exp( (dist - R) / (2. * T) ) );
  }

  // Reorders the nodes by decreasing degree. edges[0] is the node with the
  // largest degree.
  // new_permutation will contain the reordering of the nodes. permutation[i]=j
  // means that old node i is now node j.
  void sortByDegrees(vector<int>* new_permutation) override;

  // Returns the log-likelihood of a node being at the position where it
  // currently is, based on the edges to other nodes.
  // Only considers edges to nodes that lie above the given layer.
  virtual double nodeLogLikelihood(int node, int layer);

  // Computes the energy of an embedding. The smaller the energy, the better.
  // Only takes nodes into account with radial coordinate < max_rad.
  virtual double computeEnergy(int min_layer);

  // DEPRECATED: does not initialize class fields properly
  Hyperbolic(int n) : Hyperbolic(n, 1, 1, 1) {};

  // Radius R, scaling alpha, proability decay T
  Hyperbolic(int n, double R, double alpha, double T)
    : Graph(n), R(R), alpha(alpha), T(T) {};
  
  // The caller has ownership of the created objects.
  Hyperbolic subgraph(vector<bool> vertices);
  Hyperbolic giantSubgraph();

  virtual ~Hyperbolic() {};

  // Saves the edges into filename-links.txt and the hyperbolic coordinates in
  // filename-embedded.txt.
  void printToFile(const char* filename,
                   const vector<bool>* ignore_nodes = NULL,
                   bool use_labels = true) override;
};

class HyperbolicTrivial : public Hyperbolic {
 public:
  HyperbolicTrivial(int n, double R, double alpha, double T = 0.);
  
  HyperbolicTrivial(int n, vector<HYPT> pts, double R, double alpha,
                    double T = 0.);
  
  ~HyperbolicTrivial() {};
};
