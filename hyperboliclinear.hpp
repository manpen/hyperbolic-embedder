#pragma once

#include <forward_list>
#include <unordered_set>

#include "hyperbolic.hpp"
#include "geometricds.hpp"

using std::unordered_set;

class HyperbolicLinear : public Hyperbolic {
 public:
  HyperbolicLinear(int n, double R, double alpha, double T) :
    Hyperbolic(n, R, alpha, T), L(-1) {};

  static HyperbolicLinear linearSampling(int n, double R, double alpha,
                                         double T);

  // Samples edges in linear time.
  void sampleEdges();

  // Embeds the graph in the hyperbolic plane.
  void hyperbolicEmbedding();

  // See above. Only optimizes between the given layers (inclusively). Useful
  // for random restarts. If close is set, the nodes in max_layer will initially
  // be placed close to their neighbors in max_layer+1; otherwise they will be
  // distributed randomly.
  void findGoodEmbedding(int max_layer, int min_layer, bool close);

  // Computes the log-likelihood for node to be where it currently is.
  // Only consideres edges (and non-edges) to nodes in layer >= min_layer.
  double nodeLogLikelihood(int node, int min_layer) override;

  // Returns the layer of a node at radial coordinate rad.
  inline int getLayer(int node) const {
    return (int) ((R - pts[node].r) / log(4));
  }

  // Returns true if node has an active neighbor that lies in the given layer
  // or higher.
  bool hasNeighborAboveLayer(int node, int layer, const vector<bool>& active);

  HyperbolicLinear subgraph(vector<bool> vertices);
  HyperbolicLinear giantSubgraph();

  ~HyperbolicLinear() {};

  vector<unordered_set<int>> weight_layers;
  vector<PointsFromCells> points_from_cells;
  GeometricDS<TreePoints> geometric_ds;
  int L;

  long long num_edges {0};

 private:
  // Builds weight layers using the radial coordinates.
  void constructWeightLayers();

  // Builds the PointsFromCells data structure for all weight layers above
  // min_layer.
  void constructPointsFromCells(int min_layer);

  // For a given node, estimates the best angle by looking at active neighbors
  // of the node that lie in layer >= min_layer.
  // Automatically sets the angle.
  void estimateBestAngle(int node, int min_layer, const vector<bool>& active);

  // Call once before embedding to initialize core of the graph.
  // Returns the layer from where the embedding algo has to be run.
  int prepareEmbedding();
};
