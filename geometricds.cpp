//
//  geometricds.cpp
//  hyperbolic-embedder
//
//  Created by Anton Krohmer on 27.04.16.
//
//

#include <glog/logging.h>

#include "geometricds.hpp"

#define FOR(i,n) for(int (i)=0;(i)<(n);(i)++)
#define FORB(i,a,n) for(int (i)=(a);(i)<(n);(i)++)
#define FOREACH(it,c) for(__typeof((c).begin()) it=(c).begin();it!=(c).end();++it)
#define PB emplace_back
#define MP make_pair

namespace {
  inline int mod(int a, int b) {
    return ((a % b) + b) % b;
  }
}

int TreePoints::cellIntersect(int cell, int layer) {
  CHECK_GE(layer, this->layer);
  CHECK_LT(cell, tree_points[layer - this->layer].size());

  return (int)tree_points[layer - this->layer][cell].size();
}

const std::unordered_set<int>& TreePoints::cellNodes(int cell, int layer) {
  CHECK_GE(layer, this->layer);
  CHECK_LT(cell, tree_points[layer - this->layer].size());

  return tree_points[layer - this->layer][cell];
}

void TreePoints::eraseNode(int node, double phi) {
  int index = (int) (num_cells * phi / (2 * M_PI));

  FORB(i,this->layer,L) {
    CHECK(tree_points[i - this->layer][index].find(node) !=
          tree_points[i - this->layer][index].end());
    tree_points[i - this->layer][index].erase(node);

    index /= 2;
  }
}

void TreePoints::addNode(int node, double phi) {
  int index = (int) (num_cells * phi / (2 * M_PI));
  FORB(i,this->layer,L) {
    tree_points[i - this->layer][index].insert(node);
    index /= 2;
  }
}

TreePoints::TreePoints(int layer, int L, const unordered_set<int>* P,
                       const Hyperbolic& H) :
layer(layer), L(L) {
  num_cells = 1 << (L-1 - layer);

  tree_points.resize(L-layer);
  int current_cells = num_cells;
  FORB(i,layer,L) {
    tree_points[i - layer].resize(current_cells);
    current_cells /= 2;
  }

  for (int node : (*P)) {
    CHECK_LT(node, H.n);
    int index = (int) (num_cells * H.pts[node].phi / (2 * M_PI));
    FORB(i,layer,L) {
      CHECK_LT(index, tree_points[i - layer].size());
      CHECK_LT(i-layer, tree_points.size());
      tree_points[i - layer][index].insert(node);
      index /= 2;
    }
  }
}

PointsFromCells::PointsFromCells(int layer, int L, const unordered_set<int>* P,
                                 const Hyperbolic& H)
:  P(P), layer(layer), L(L) {

  FOR(i,H.n)
  CHECK(H.pts[i].phi >= 0) << i << ": " << H.pts[i].phi;

  num_cells = 1 << (L-1 - layer);
  //  LOG(INFO) << "layer: " << layer << ", num cells: " << num_cells;

  num_points_in_cell.resize(num_cells);
  for (int i : (*P)) {
    num_points_in_cell[(int) (num_cells * (H.pts[i].phi / (2*M_PI)))] += 1;
    //    LOG(INFO) << "Added point at " << H.pts[i].phi << " to cell "
    //        << (int) (num_cells * (H.pts[i].phi / (2*M_PI)))
    //        << " where num_cells = " << num_cells;
  }

  s.resize(num_cells);
  s[0] = 0;
  FORB(i,1,num_cells)
  s[i] = s[i-1] + num_points_in_cell[i-1];

  // A.clear();
  A.resize(P->size(),-1);
  vector<int> pointers(num_cells, 0);
  for (int i : (*P)) {
    int cell = (int) (num_cells * (H.pts[i].phi / (2*M_PI)));
    CHECK_EQ(A[s[cell] + pointers[cell]], -1);
    A[s[cell] + pointers[cell]] = i;
    pointers[cell] += 1;
  }
}

int PointsFromCells::cellIntersect(int cell, int layer) const {
  CHECK_GE(layer, this->layer);
  CHECK_GE(cell, 0);
  int index_begin = cell << (layer - this->layer);
  int index_end = ((cell+1) << (layer - this->layer)) - 1;
  //  LOG(INFO) << "This layer: " << this->layer << ", other layer: " << layer
  //            << ", cell: " << cell << " (index " << index_begin << "-" << index_end << ")";

  CHECK_LT(index_end, num_cells);
  return s[index_end] - s[index_begin] + num_points_in_cell[index_end];
}

int PointsFromCells::kthPoint(int cell, int layer, int k) const {
  CHECK_GE(layer, this->layer);
  int index_begin = cell << (layer - this->layer);
  CHECK_LT(k,cellIntersect(cell, layer));
  CHECK_LT(s[index_begin] + k, A.size());
  CHECK_GE(A[s[index_begin] + k], 0);
  return A[s[index_begin] + k];
}

LayerPartitioning::LayerPartitioning(int i, int j, int L) :
L(L), has_next(true), next_type(true), current_cell(0), current_pair_cell(0) {
  current_layer = i+j;
  if (current_layer >= L-1)
    current_layer = L-1;
}

//LayerPartitioning::LayerPartitioning(int i, int node, int L,
//                                     const HyperbolicLinear& H) :
//L(L), has_next(true), next_type(true) {
//  current_layer = i + (int) ((H.R - H.pts[node].r) / log(4));
//  if (current_layer > L-1)
//    current_layer = L-1;
//
//  current_cell = (int) ((1 << (L-1 - current_layer)) *
//                        H.pts[node].phi / (2 * M_PI));
//  current_pair_cell = current_cell;
//}

LayerPartitioning::LayerPartitioning(int i, int node, int L,
                                     const Hyperbolic& H, double R)
: L(L), has_next(true), next_type(true) {
  current_layer = i + (int)((R - H.pts[node].r) / log(4));
  if (current_layer > L - 1) current_layer = L - 1;

  current_cell =
  (int)((1 << (L - 1 - current_layer)) * H.pts[node].phi / (2 * M_PI));
  current_pair_cell = current_cell;
}

void LayerPartitioning::getNext(int* cell, int* layer) {
  int num_cells = (1 << (L-1 - current_layer));

  *cell = mod(current_pair_cell, num_cells);
  *layer = current_layer;

  if (next_type) { // Type I
    if (num_cells == 1)
      has_next = false;
    else if (current_pair_cell == current_cell)
      current_pair_cell = current_cell - 1;
    else if (current_pair_cell == current_cell - 1 && num_cells > 2)
      current_pair_cell = current_cell + 1;
    else if (num_cells > 2) {
      next_type = false;
      current_pair_cell = current_cell + 2;
    } else
      has_next = false;
  } else { // Type II
    CHECK_GE(num_cells, 4);
    if (num_cells == 4)
      has_next = false;
    else if (current_pair_cell == current_cell - 2) {
      current_layer += 1;
      current_cell /= 2;
      current_pair_cell = current_cell + 2;
    } else if (current_cell % 2 == 0) {
      if (current_pair_cell == current_cell + 2)
        current_pair_cell = current_cell + 3;
      else if (current_pair_cell == current_cell + 3)
        current_pair_cell = current_cell - 2;
      else
        CHECK(false);
    } else {
      if (current_pair_cell == current_cell + 2)
        current_pair_cell = current_cell - 3;
      else if (current_pair_cell == current_cell - 3)
        current_pair_cell = current_cell - 2;
      else
        CHECK(false);
    }
  }
}


void LayerPartitioning::getNext(int* cellA, int* layerA,
                                int* cellB, int* layerB) {
  int num_cells = (1 << (L - 1 - current_layer));

  *cellA  = current_cell;
  *cellB = mod(current_pair_cell,num_cells);
  *layerA = current_layer;
  *layerB = current_layer;

  CHECK_LT(*cellA, num_cells);
  CHECK_LT(*cellB, num_cells);

  if (next_type) { // Type I
    if (num_cells == 1)
      has_next = false;
    else if (current_pair_cell == current_cell)
      current_pair_cell = current_cell - 1;
    else if (current_pair_cell == current_cell - 1 && num_cells > 2)
      current_pair_cell = current_cell + 1;
    else if (current_cell < num_cells - 1)
      ++current_cell;
    else if (num_cells > 2) {
      next_type = false;
      current_cell = 0;
      current_pair_cell = 2;
    } else
      has_next = false;
  } else { // Type II
    if (current_cell % 2 == 0) {
      if (num_cells == 4) {  // if num cells == 4, the case is degenerate
        current_cell = 1;
        current_pair_cell = 3;
      } else if (current_pair_cell == current_cell + 2)
        current_pair_cell = current_cell + 3;
      else if (current_pair_cell == current_cell + 3)
        current_pair_cell = current_cell - 2;
      else
        current_cell = current_cell + 1;
    } else {
      if (num_cells == 4)
        has_next = false;
      else if (current_pair_cell == current_cell - 3)
        current_pair_cell = current_cell - 2;
      else if (current_pair_cell == current_cell - 2)
        current_pair_cell = current_cell + 2;
      else if (current_cell < num_cells - 1) {
        current_cell += 1;
        current_pair_cell = current_cell + 2;
      } else if (current_layer < L-2) { // L-1 contains one cell, ergo L-2
        // cells are adjacent
        current_layer += 1;
        current_cell = 0;
        current_pair_cell = 2;
      } else {
        has_next = false;
      }
    }
  }
}


template<typename T>
void GeometricDS<T>::init(Hyperbolic* H) {
  m_G = H;
  m_R = H->R;
  m_L = ((int)(m_R / log(4)) + 1);
  // separate vertices into L layers
  vector<unordered_set<int>> weight_layers(m_L);
  for (int v = 0; v < m_G->n; ++v) weight_layers[getLayer(v)].insert(v);

  // create the tree nodes
  m_tree_points.clear();
  for (int i = 0; i < m_L; ++i)
    m_tree_points.emplace_back(T(i, m_L, &weight_layers[i], *m_G));
}

template<typename T>
void GeometricDS<T>::forAllNodesDo(
    int u, function<void(int v, int count, HYPT pos, bool type)> f,
    int min_layer, bool long_distances) {
  int cell, layer;
  for (int i = min_layer; i < m_L; ++i) {
    LayerPartitioning lp(i, u, m_L, *m_G, m_R);
    while (lp.hasNext()) {
      if (lp.nextType()) {  // vertex close to v1
        lp.getNext(&cell, &layer);

        applyToPointsInCell(u, cell, layer, m_tree_points[i], f);

      } else {  // vertex far away from v1
        if (!long_distances) return;
        lp.getNext(&cell, &layer);

        // number of points in this layer
        int cell_points = m_tree_points[i].cellIntersect(cell, layer);
        if (cell_points == 0) continue;

        // the average point in this layer
        HYPT avg_pt;
        avg_pt.phi = 2 * M_PI * (cell + 0.5) / (1 << (m_L - 1 - layer));
        avg_pt.r = m_R - (i + 0.5) * log(4);

        // compute force as if v1 and all points in the current cell
        // are not connected
        f(-1, cell_points, avg_pt, false);
      }
    }
  }
}

template <typename T>
void GeometricDS<T>::applyToPointsInCell(
    int u, int cell, int layer, TreePoints& tp,
    function<void(int v, int count, HYPT pos, bool type)> f) {
  for (int v : tp.cellNodes(cell, layer)) {
    if (u != v) f(v, 1, m_G->pts[v], true);
  }
}

template <typename T>
void GeometricDS<T>::applyToPointsInCell(
    int u, int cell, int layer, PointsFromCells& pfc,
    function<void(int v, int count, HYPT pos, bool type)> f) {
  int max = pfc.cellIntersect(cell, layer);
  for (int i = 0; i < max; ++i) {
    int v = pfc.kthPoint(cell, layer, i);
    if (u != v) f(v, 1, m_G->pts[v], true);
  }
}

template<>
void GeometricDS<TreePoints>::move(int v, HYPT new_pos) {
  // remove v from the data structure
  m_tree_points[getLayer(v)].eraseNode(v, m_G->pts[v].phi);

  // move v to the new position
  m_G->pts[v].r = new_pos.r;
  m_G->pts[v].phi = new_pos.phi;

  // reinsert v
  m_tree_points[getLayer(v)].addNode(v, m_G->pts[v].phi);
}

template<typename T>
void GeometricDS<T>::move(int v, HYPT new_pos) {
  // remove v from the data structure
  // m_tree_points[getLayer(v)].eraseNode(v, m_G.pts[v].phi);
  
  // move v to the new position
  m_G->pts[v].r = new_pos.r;
  m_G->pts[v].phi = new_pos.phi;
  
  // reinsert v
  // m_tree_points[getLayer(v)].addNode(v, new_pos.phi);
}

template class GeometricDS<TreePoints>;
template class GeometricDS<PointsFromCells>;