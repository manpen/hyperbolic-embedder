#include <cmath>
#include <cstdio>
#include <fstream>
#include <glog/logging.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_monte_plain.h>
#include <gsl/gsl_monte_miser.h>
#include <gsl/gsl_monte_vegas.h>
#include <iostream>
#include <list>
#include <limits>
#include <sstream>
#include <vector>
#include <algorithm>
#include <unordered_map>

#include "graph.hpp"
#include "random.hpp"
#include "hyperbolic.hpp"
#include "hyperboliclinear.hpp"
#include "NLEHelper.hpp"
#include "springembedder.hpp"

#define FOR(i,n) for(int (i)=0;(i)<(n);(i)++)
#define FORB(i,a,n) for(int (i)=(a);(i)<(n);(i)++)
#define FOREACH(it,c) for(__typeof((c).begin()) it=(c).begin();it!=(c).end();++it)
#define PB emplace_back
#define MP make_pair

using std::vector;

namespace {
  inline double mod(double a, double b) {
    return fmod(fmod(a, b) + b, b);
  }
}

Hyperbolic* HyperbolicLinear::subgraph(vector<bool> vertices) {
  vector<int> index (n+1);
  index[0] = 0;
  FORB(i,1,n+1) index[i] = index[i-1] + (vertices[i-1] ? 1 : 0);
  int m = index[n];
  HyperbolicLinear* H = new HyperbolicLinear(m, R, alpha, T);
  FOR(i,n)
  if (vertices[i]) {
    H->labels[index[i]] = labels[i];
    H->pts.PB(pts[i]);
    FOR(j,(int)edges[i].size()) {
      int jj = edges[i][j];
      if (vertices[jj])
        H->edges[index[i]].PB(index[jj]);
    }
  }

  return H;
}

Hyperbolic* HyperbolicLinear::giantSubgraph() {
  vector<int> comp = components();
  vector<bool> gcc = giant(comp);
  return subgraph(gcc);
}

HyperbolicLinear* HyperbolicLinear::linearSampling(int n, double R,
                                                   double alpha, double T) {
  HyperbolicLinear* H = new HyperbolicLinear(n, R, alpha, T);
  H->pts.resize(n);
  FOR(i,n)
    H->pts[i] = HYPT::randomHYPT(R, alpha);
  H->sampleEdges();

  return H;
}

void HyperbolicLinear::constructWeightLayers() {
  weight_layers.clear();
  // determine weight layers
  L = (int) (R / log(4)) + 1;
  CHECK_LE(L,std::numeric_limits<int>::digits);
  weight_layers.resize(L);
  FOR(i,n) {
    weight_layers[getLayer(i)].insert(i);
  }
}

void HyperbolicLinear::constructPointsFromCells(int min_layer) {
  points_from_cells.clear();
  points_from_cells.resize(min_layer);
  FORB(i,min_layer,L)
  points_from_cells.PB(PointsFromCells(i, L, &weight_layers[i], *this));
}

void HyperbolicLinear::sampleEdges() {
  constructWeightLayers();
  constructPointsFromCells(0);

  // NOTE: We could do this using the geometric DS, but this is
  // a factor log n faster

  int cellA, layerA, cellB, layerB;

  FOR(i,L) {
    FORB(j,i,L) {
      // construct partitioning

      LayerPartitioning lp(i,j,L);
      while (lp.hasNext()) {
        // iterate over all pairs of cells that are in layers > i+j

        if (lp.nextType()) {
          // type I
          lp.getNext(&cellA, &layerA, &cellB, &layerB);

          int cellA_nodes = points_from_cells[i].cellIntersect(cellA, layerA);
          int cellB_nodes = points_from_cells[j].cellIntersect(cellB, layerB);
          FOR(x,cellA_nodes) {
            FOR(y,cellB_nodes) {
              int u = points_from_cells[i].kthPoint(cellA, layerA, x);
              int v = points_from_cells[j].kthPoint(cellB, layerB, y);

              if ((j > i || u > v) &&
                  randdbl() < Hyperbolic::prob(HYPT::dist(pts[u], pts[v]),
                                               R, T)) {
                    edges[u].PB(v);
                    edges[v].PB(u);
                  }
            }
          }
        } else {
          // type II
          if (T == 0)
            break; // Step model, in that case there are no edges of type II

          lp.getNext(&cellA, &layerA, &cellB, &layerB);

          int diff = std::abs((1 << (L - 2 - layerA)) -
                              std::abs((1 << (L - 2 - layerA)) -
                                       std::abs(cellA - cellB))) - 1;
          CHECK_GT(diff, 0);
          CHECK_LT(diff, 3);

          double min_phi_dist = 2 * M_PI * diff / (1 << (L-1-layerA));
          double rad1 = R - (i+1) * log(4);
          double rad2 = R - (j+1) * log(4);
          double p_bar = Hyperbolic::prob(HYPT::dist(rad1, rad2, 0,
                                                     min_phi_dist), R, T);

          long r=0;
          if (p_bar == 0)
            continue;
          else if (p_bar == 1)
            r += 1;
          else {
            double rd = (std::ceil(log(randdbl()) / log(1-p_bar))); // r >= 1
            if (!std::isfinite(rd) || rd > std::numeric_limits<long>::max())
              continue;
            r += (long) rd;
          }

          CHECK_GT(r,0);
          int cellA_points = points_from_cells[i].cellIntersect(cellA,layerA);
          int cellB_points = points_from_cells[j].cellIntersect(cellB,layerB);
          if (1 <= r && r <= cellA_points * cellB_points) {
            // determine r-th pair in V_i^A x V_j^B

            int u = points_from_cells[i].kthPoint(cellA, layerA,
                                                  (r-1) % cellA_points);
            int v = points_from_cells[j].kthPoint(cellB, layerB,
                                                  (int)((r-1) / cellA_points));

            double p = Hyperbolic::prob(HYPT::dist(pts[u],
                                                   pts[v]),
                                        R, T);
            CHECK_LE(p,p_bar) << "deltaphi: " << pts[u].phi - pts[v].phi
            << ", r_u: " << pts[u].r << ", r_v: " << pts[v].r;

            if ((j > i || u > v) && randdbl() < p / p_bar) {
              edges[u].PB(v);
              edges[v].PB(u);
            }

            if (p_bar == 0)
              continue;
            else if (p_bar == 1)
              r += 1;
            else {
              double rd = (std::ceil(log(randdbl()) / log(1-p_bar))); // r >= 1
              if (!std::isfinite(rd) || rd > std::numeric_limits<long>::max())
                continue;
              r += (long) rd;
            }
          }

        }
      }
    }
  }
}


double HyperbolicLinear::nodeLogLikelihood(int node, int min_layer) {

  // enables constant time lookup for neighbors
  std::unordered_set<int> neighbors;

  double loglikelihood = 0;
  for (int neighbor : edges[node]) {
    neighbors.insert(neighbor);

    if (getLayer(neighbor) >= min_layer) {
      double dist = HYPT::dist(pts[neighbor], pts[node]);
      double prob = Hyperbolic::prob(dist, R, T);
      if (prob > 0)
        loglikelihood += log(prob);
      else
        loglikelihood += - (dist - R)/(2 * T);

      CHECK(std::isfinite(loglikelihood)) << neighbor
      << "-" << node << " (" << HYPT::dist(pts[neighbor], pts[node]) << ", "
      << R << ")";
    }
  }

  geometric_ds.forAllNodesDo(node, [&](int u, int count, HYPT pos, bool type) {
    if (u != node && neighbors.count(u) == 0) {
      double likelihood = log(1 - prob(HYPT::dist(pts[u], pts[node]), R, T));
      if (!std::isfinite(likelihood)) // precision too low, use estimate
        likelihood = (HYPT::dist(pts[u], pts[node]) - R) / (2*T);
      loglikelihood += likelihood;
    }
  }, min_layer, false);

  return loglikelihood;
}

void HyperbolicLinear::hyperbolicEmbedding() {
  CHECK_EQ(sizeOfGiant(giant(components())),n)
      << "Can only embed connected graphs.";

  vector<vector<int>> orig_edges = edges;

  constructWeightLayers();

  int max_layer = prepareEmbedding();
  geometric_ds.init(this);
  findGoodEmbedding(max_layer, 0, true);
  edges = orig_edges;
}

void HyperbolicLinear::estimateBestAngle(int node, int min_layer,
                                         const vector<bool>& active) {

  // compute average angle
  double sinsum = 0;
  double cossum = 0;
  for (int j : edges[node]) {
    if (getLayer(j) >= min_layer && active[j]) {
      sinsum += (1 << (L-1 - getLayer(j))) * sin(pts[j].phi);
      cossum += (1 << (L-1 - getLayer(j))) * cos(pts[j].phi);
    }
  }

  if (cossum != 0)
    pts[node].phi = mod(atan(sinsum / cossum) + (cossum < 0 ? M_PI : 0),
                        2*M_PI);
  else {
    LOG(WARNING) << "Using random angle for node: " << node;
    pts[node].phi = 2*M_PI*randdbl();
  }
}

int HyperbolicLinear::prepareEmbedding() {
  int threshold = 3;
  CHECK_GT(n, 2*threshold) << "Graph is too small to embed.";

  vector<int> timestamp(n,-1);
  int max_layer = L/2;
  int inner_nodes = 0;

  FOR(i, n) {
    if (getLayer(i) >= max_layer)
      ++inner_nodes;
    else
      break;
  }

  while (inner_nodes < 2* threshold) {
    inner_nodes = 0;
    --max_layer;
    FOR(i, n) {
      if (getLayer(i) >= max_layer)
        ++inner_nodes;
      else
        break;
    }
  }

  // Count common neighbors
  vector<vector<double>> common_neighbors(inner_nodes);
  FOR(node,inner_nodes)
  common_neighbors[node].resize(inner_nodes, 0);

  FOR(node, inner_nodes) {
    for (int neighbor : edges[node])
      timestamp[neighbor] = node;

    FORB (othernode, node+1, inner_nodes)
    for (int neighbor : edges[othernode])
      if (timestamp[neighbor] == node)
        ++common_neighbors[node][othernode];
  }


  // compute the angle between two nodes based on the number of common neighbors
  FOR(i, inner_nodes) {
    FORB(j, i+1, inner_nodes) {
      int smalldeg = (int) edges[i].size();
      int largedeg = (int) edges[j].size();
      if (smalldeg > largedeg) std::swap(smalldeg, largedeg);
      common_neighbors[i][j] = largedeg *
      pow(smalldeg / common_neighbors[i][j], 1/(2 * alpha - 1)) / n;
      common_neighbors[j][i] = common_neighbors[i][j];
    }
  }

  // For each inner node, keep only the edges (to other inner nodes) with
  // the highest score

  vector<vector<int>> important_edges(inner_nodes);
  FOR(i, inner_nodes)
    FOR(j, inner_nodes)
      if (std::count(edges[i].begin(), edges[i].end(), j) > 0)
        important_edges[i].PB(j);

  FOR(i, inner_nodes) {
    std::sort(important_edges[i].begin(), important_edges[i].end(),
              [&](int n1, int n2) {
                return common_neighbors[i][n1] < common_neighbors[i][n2];
              });
  }

  FOR(i, inner_nodes)
    if (threshold < important_edges[i].size())
      important_edges[i].resize(threshold);
  FOR(i, inner_nodes)
    for (int inner_neighbor : important_edges[i])
      if (std::count(important_edges[inner_neighbor].begin(),
                     important_edges[inner_neighbor].end(), i) == 0)
        important_edges[inner_neighbor].PB(i);

  FOR(node, inner_nodes) {
    vector<int> new_neighbors;
    // remove edges
    for (int inner_neighbor : important_edges[node])
      new_neighbors.PB(inner_neighbor);

    for (int neighbor : edges[node])
      if (getLayer(neighbor) < max_layer) new_neighbors.PB(neighbor);

    edges[node].swap(new_neighbors);
  }


  vector<HYPT> orig_pts(pts);
  SpringEmbedder spring(*this);

  // Restart a few times to get a good embedding
  vector<HYPT> best_pts(pts);
  double best_score = std::numeric_limits<double>::max();
  FOR(i,5) {
    spring.run(common_neighbors, 600, 1, 0, true);
    if (spring.score(common_neighbors) < best_score) {
      best_score = spring.score(common_neighbors);
      best_pts = pts;
    }
  }

  FOR (i, n) {
    pts[i].r = orig_pts[i].r;
    pts[i].phi = best_pts[i].phi;
    if (std::isnan(pts[i].phi))
      pts[i].phi = 0;
  }

  return max_layer - 1;
}

bool HyperbolicLinear::hasNeighborAboveLayer(int node, int layer,
                                             const vector<bool>& active) {
  for (int neighbor : edges[node])
    if (getLayer(neighbor) >= layer && active[neighbor])
      return true;

  return false;
}

void HyperbolicLinear::findGoodEmbedding(int max_layer, int min_layer,
                                         bool close) {
  CHECK_GE(max_layer, min_layer);
  double lastEnergy = std::numeric_limits<double>::lowest();
  double improvement = 1;
  vector<bool> active(n, false);
  vector<int> smallest_neighbor(n,-1);
  FOR(i,n) {
    if (getLayer(i) > max_layer || (!close && getLayer(i) == max_layer)) {
      active[i] = true;
      for (int neighbor : edges[i])
        if (smallest_neighbor[neighbor] < 0 ||
            pts[i].r > pts[smallest_neighbor[neighbor]].r)
          smallest_neighbor[neighbor] = i;
    }
  }


  // proceed iteratively in L/2 = log n steps
  // take all nodes in layer >= L/2 - i
  for (int layer = max_layer; layer >= min_layer; --layer) {
    improvement = 1;
    lastEnergy = std::numeric_limits<double>::lowest();

    // Place new arriving nodes close to their neighbors in the previous layers.
    if (layer < max_layer || close) {
      FOR(i,n) {
        if (getLayer(i) == layer && hasNeighborAboveLayer(i, layer+1, active)) {
          double old_phi = pts[i].phi;
          estimateBestAngle(i, layer+1, active);
          double new_phi = pts[i].phi;
          pts[i].phi = old_phi;
          geometric_ds.move(i, HYPT(pts[i].r, new_phi));
          active[i] = true;
          for (int neighbor : edges[i])
            if (smallest_neighbor[neighbor] < 0 ||
                pts[i].r > pts[smallest_neighbor[neighbor]].r)
              smallest_neighbor[neighbor] = i;
        }
      }
    }

    while(improvement > 0.001) {
      FOR(i,n) {
        if (getLayer(i) >= layer && !active[i]
            && hasNeighborAboveLayer(i, layer, active)) {
          double old_phi = pts[i].phi;
          estimateBestAngle(i, layer, active);
          double new_phi = pts[i].phi;
          pts[i].phi = old_phi;
          geometric_ds.move(i, HYPT(pts[i].r, new_phi));
          active[i] = true;
          for (int neighbor : edges[i])
            if (smallest_neighbor[neighbor] < 0 ||
                pts[i].r > pts[smallest_neighbor[neighbor]].r)
              smallest_neighbor[neighbor] = i;
        }

        if (!active[i]) continue;

        double old_phi = pts[i].phi;
        double best_phi = pts[i].phi;
        double best_ll = nodeLogLikelihood(i, layer);

        estimateBestAngle(i, layer, active);
        double best_angle_est = pts[i].phi;
        if (best_ll < nodeLogLikelihood(i, layer)) {
          best_phi = pts[i].phi;
          best_ll = nodeLogLikelihood(i, layer);
        }

        // retry a few times to find the best spot
        FOR(x,10 + log(n)) {
          pts[i].phi = mod(
              best_angle_est + randgaussian(
                   4 * exp((R-pts[i].r-pts[smallest_neighbor[i]].r)/2.0)),
                           2*M_PI);

          if (best_ll < nodeLogLikelihood(i, layer)) {
            best_phi = pts[i].phi;
            best_ll = nodeLogLikelihood(i, layer);
          }
        }

        pts[i].phi = old_phi;
        geometric_ds.move(i, HYPT(pts[i].r, best_phi));
      }
      
      // compute total energy
      double energy = computeEnergy(layer);
      improvement = (lastEnergy-energy)/lastEnergy;
      
      lastEnergy = energy;
      LOG(INFO) << "Layer " << layer << " has log-likelihood " << lastEnergy
      << " (Improvement by " << improvement << ")";
    }
  }
}
