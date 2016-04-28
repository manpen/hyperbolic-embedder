#include <vector>
#include <list>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <glog/logging.h>
#include <set>
#include <queue>
#include <algorithm>
#include <unordered_map>

#include "graph.hpp"
#include "random.hpp"
#include "hyperbolic.hpp"

#define FOR(i,n) for(int (i)=0;(i)<(n);(i)++)
#define FORB(i,a,n) for(int (i)=(a);(i)<(n);(i)++)
#define FOREACH(it,c) for(__typeof((c).begin()) it=(c).begin();it!=(c).end();++it)
#define PB push_back
#define MP make_pair

using std::pair;
using std::max;
using std::min;
using std::queue;
using std::vector;
using std::set;
using std::make_pair;
using std::string;
using std::isfinite;

HYPT::HYPT(double r, double phi) {
  this->r = r;
  this->phi = phi;
}

HYPT HYPT::randomHYPT(double R, double alpha) {
  double r = acosh( randdblpos() * (cosh(alpha * R) - 1.) + 1. ) / alpha;
  return HYPT(r, randdblpos() * 2. * M_PI);
}

double HYPT::dist(HYPT p, HYPT q) {
  return dist(p.r, q.r, p.phi, q.phi);
}

double HYPT::dist(double r1, double r2, double phi1, double phi2) {
  return acosh( max( 1., cosh(r1 - r2) + (1. - cos(phi1 - phi2))*sinh(r1)*sinh(r2) ) );
}

double Hyperbolic::nodeLogLikelihood(int node, int layer) {
  double ll = 0;
  double rad_max = R - layer * log(4);
  FOR(i,n) {
    if (i != node && pts[i].r <= rad_max) {
      double dist = HYPT::dist(pts[i], pts[node]);
      double p = prob(dist, R, T);
      if (adjacent(i,node)) {
        if (p > 0)
          ll += log(p);
        else
          ll += - (dist - R)/(2 * T);

        CHECK(isfinite(ll)) << i
        << "-" << node << " (" << HYPT::dist(pts[i], pts[node]) << ", "
        << R << ")";
      } else {
        double likelihood = log(1 - p);
        if (!isfinite(likelihood)) // precision too low, use estimate
          likelihood = (dist - R) / (2*T);
        ll += likelihood;
      }
    }
  }
  return ll;
}

double Hyperbolic::computeEnergy(int min_layer) {
  double energy = 0;
  double max_rad = R - min_layer * log(4);
  FOR(i,n) {
    if (pts[i].r > max_rad) continue;
    energy += nodeLogLikelihood(i,min_layer);
  }

  return energy;
}

HyperbolicTrivial::HyperbolicTrivial(int n, vector<HYPT> pts, double R,
                                     double alpha, double T)
: Hyperbolic(n, R, alpha, T) {
  this->n = n;
  this->R = R;
  this->T = T;
  this->pts = pts;
  FOR(i,n) {
    FOR(j,i) {
      if (randdblpos() <= prob(HYPT::dist(pts[i],pts[j]), R, T)) {
        edges[j].PB(i);
        edges[i].PB(j);
      }
    }
  }
}

Hyperbolic* Hyperbolic::subgraph(vector<bool> vertices) {
  vector<int> index (n+1);
  index[0] = 0;
  FORB(i,1,n+1) index[i] = index[i-1] + (vertices[i-1] ? 1 : 0);
  int m = index[n];
  Hyperbolic* G = new Hyperbolic(m, R, alpha, T);
  FOR(i,n)
  if (vertices[i]) {
    G->labels[index[i]] = labels[i];
    G->pts.PB(pts[i]);
    FOR(j,(int)edges[i].size()) {
      int jj = edges[i][j];
      if (vertices[jj])
        G->edges[index[i]].PB(index[jj]);
    }
  }
  return G;
}

Hyperbolic* Hyperbolic::giantSubgraph() {
  vector<int> comp = components();
  vector<bool> gcc = giant(comp);
  return subgraph(gcc);
}

void Hyperbolic::printToFile(const char *filename,
                             const vector<bool>* ignore_nodes,
                             bool use_labels) {
  std::string filename_appended(filename);
  filename_appended += "-links.txt";

  Graph::printToFile(filename_appended.c_str(), ignore_nodes);

  filename_appended.assign(filename);
  filename_appended += "-coordinates.txt";

  std::ofstream file;
  file.open(filename_appended.c_str());

  file << "n\tR\talpha\tT" << std::endl;
  file << n << "\t" << R << "\t" << alpha << "\t" << T << std::endl;

  //file.flags(std::ios::scientific);
  file.precision(std::numeric_limits<double>::digits10 + 1);

  FOR(i,n) {
    if (ignore_nodes == NULL || !ignore_nodes->at(i)) {
      if (use_labels) file << labels[i];
      else            file << i;
      file << "\t" << pts[i].r << "\t" << pts[i].phi * 360 / (2*M_PI)
      << std::endl;
    }
  }

  file.close();
}

void Hyperbolic::sortByDegrees(vector<int>* new_permutation) {
  Graph::sortByDegrees(new_permutation);

  vector<HYPT> new_pts(n);
  FOR(i, n) {
    new_pts[new_permutation->at(i)] = pts[i];
  }

  pts.swap(new_pts);
}
