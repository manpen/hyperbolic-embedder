//
//  geometricds.hpp
//  hyperbolic-embedder
//
//  Created by Anton Krohmer on 27.04.16.
//
//

#pragma once

#include <functional>
#include <unordered_set>

#include "hyperbolic.hpp"

using std::unordered_set;

using std::function;

class PointsFromCells {
public:
  // Dummy constructor for speedup
  PointsFromCells() : P(NULL), layer(-1), L(-1), num_cells(-1) {}

  PointsFromCells(int layer, int L, const unordered_set<int>* P,
                  const Hyperbolic& H);

  // Returns the number of points that lie in P and the given cell.
  // The layer has 2^(layer) many cells.
  int cellIntersect(int cell, int layer) const ;

  // Returns the k-th point in C \cap P (in a fixed geometric ordering of
  // C \cap P depending only on P and volume).
  int kthPoint(int cell, int layer, int k) const;

private:
  // Stores the point set.
  const unordered_set<int>* P;

  // The layer from which we allow cells. Each layer has 2^(L-layer) many
  // cells.
  int layer;
  int L;
  int num_cells;

  // Stores the number of points in each cell.
  vector<int> num_points_in_cell;

  // Prefix sum where s[i] containing the number of points in P that appear in
  // cells j < i
  vector<int> s;

  // A[s[i]+k] contains the k-th point in cell i.
  vector<int> A;
};

class TreePoints {
public:
  TreePoints(int layer, int L, const unordered_set<int>* P,
             const Hyperbolic& H);

  int cellIntersect(int cell, int layer);

  const std::unordered_set<int>& cellNodes(int cell, int layer);

  // Moves a node to another cell.
  void eraseNode(int node, double phi);
  void   addNode(int node, double phi);

private:
  int layer;
  int L;
  int num_cells;

  // TODO could save space by only storing each node once and constructing an
  // iterator that goes over all child cells.
  vector<vector<std::unordered_set<int>>> tree_points;
};

class LayerPartitioning {
public:
  // The two layers i,j for which we partition
  LayerPartitioning(int i, int j, int L);

//  // Gives you a list of cell pairs when one node is fixed
//  LayerPartitioning(int i, int node, int L, const HyperbolicLinear& H);

  LayerPartitioning(int i, int node, int L, const Hyperbolic& G, double R);

  // Returns true if there are unchecked pairs of cells
  bool hasNext() { return has_next; }

  // Returns true for type I and false for type II
  bool nextType() { return next_type; }

  // Returns the next pair of cells
  void getNext(int* cellA, int* layerA, int* cellB, int* layerB);

  // Returns the next cell (if one node is fixed)
  void getNext(int* cell, int* layer);

private:
  int L;
  bool has_next;
  bool next_type;
  int current_cell;
  int current_pair_cell;
  int current_layer;
};

template<typename T>
class GeometricDS {
public:
  // compute the data structure with the current positions
  // and radius R
  void init(Hyperbolic* G);

  // Iterate over all nodes v from the perspective of u and apply the
  // function f.  There are two cases.
  // 1. v is a vertex close to u, then count is 1, pos is the current
  //    position of v and type is true
  // 2. v represents a set of vertices, then v is -1, count is the
  //    number of vertices v represents, pos is an estimated position
  //    and type is false. This case is only handled if long_distances=true.
  void forAllNodesDo(int u,
                     function<void(int v, int count, HYPT pos, bool type)> f,
                     int min_layer=0, bool long_distances=true);

  // moves the vertex v to new_pos while keeping the data structure up
  // to date
  void move(int v, HYPT new_pos);

private:
  // get layer for a given radius
  inline int getLayer(double r) { return (int)((m_R - r) / log(4)); }

  // get layer for a given vertex
  inline int getLayer(int v) { return getLayer(m_G->pts[v].r); }

  void applyToPointsInCell(
                           int u, int cell, int layer, TreePoints& tp,
                           function<void(int v, int count, HYPT pos, bool type)> f);

  void applyToPointsInCell(
                           int u, int cell, int layer, PointsFromCells& pfc,
                           function<void(int v, int count, HYPT pos, bool type)> f);

  Hyperbolic* m_G;

  double m_R;
  
  // the number of layers
  int m_L;
  
  // the tree points
  vector<T> m_tree_points;
};