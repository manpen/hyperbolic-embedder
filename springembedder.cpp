
#include "springembedder.hpp"

#include <sstream>
#include <cmath>
#include <iomanip>
#include <glog/logging.h>
#include <type_traits>
#include <limits>
#include <algorithm>
#include <cstdlib>
#include <fstream>

#include "NLEHelper.hpp"
#include "random.hpp"

using std::cout;
using std::endl;
using std::isfinite;
using std::abs;
using std::max;
using std::min;
using std::floor;

SpringEmbedder::SpringEmbedder(Hyperbolic& G)
    : m_G(G),
      m_n(m_G.n),
      m_T(m_G.T),
      m_R(m_G.R),
      m_gds() {
}

void SpringEmbedder::run(int steps, int start, int pictures,
                         bool rand_init_coords) {
  // start with random positions
  if (rand_init_coords) randomCoordinates();

  // initialize the geometric data structure
  m_gds.init(&m_G);

  // initialize run
  m_steps = steps;
  m_start = start;
  m_i = start - 1;
  m_pictures = pictures;

  // run the iterations
  vector<HYPT> forces(m_n);
  for (m_i = start; m_i <= m_steps; ++m_i) {
    if (m_n == m_G.n) {
      computeForces(forces);
    } else {
      computeForcesSlow(forces);
    }
    applyForces(forces, 1.0 * m_R * cooling(), 1.0 * M_PI * cooling());
  }
}

void SpringEmbedder::run(int n, int steps, int start, int pictures,
                         bool rand_init_coords) {
  // compute the adjacency matrix for the considered subgraph
  m_adj_matrix = vector<vector<bool>>(n, vector<bool>(n, false));
  for (int u = 0; u < n; ++u) {
    for (int v : m_G.edges[u]) {
      if (v >= n) continue;
      m_adj_matrix[u][v] = true;
    }
  }

  m_n = n;
  run(steps, start, pictures, rand_init_coords);
  m_n = m_G.n;
}

void SpringEmbedder::run(const vector<vector<double>>& weights, int steps,
                         int start, int pictures, bool rand_init_coords) {
  // initialize run
  m_steps = steps;
  m_start = start;
  m_i = start - 1;
  m_pictures = pictures;
  m_n = (int) weights.size();

  if (rand_init_coords)
    for (int v = 0; v < m_n; ++v) m_G.pts[v].phi = randdblpos() * 2.0 * M_PI;

  double max_weight;
  double min_weight;
  minMaxWeights(weights, min_weight, max_weight);


  vector<double> forces(m_n);
  vector<double> abs_forces(m_n);
  for (m_i = start; m_i <= m_steps; ++m_i) {
    computeForces1D(forces, weights, min_weight, max_weight);
    applyForces1D(forces, m_max_step_1D * cooling());
  }

  m_n = m_G.n;
}

double SpringEmbedder::score(const vector<vector<double>>& weights) {
  int old_n = m_n;
  m_n = (int) weights.size();

  double max_weight;
  double min_weight;
  minMaxWeights(weights, min_weight, max_weight);

  vector<double> abs_forces(m_n);
  double score = 0.0;
  
  computeAbsForces1D(abs_forces, weights, min_weight, max_weight);

  for (int v = 0; v < m_n; ++v) {
    score += abs_forces[v];
  }

  m_n = old_n;

  return score;
}

void SpringEmbedder::randomCoordinates() {
  for (int v = 0; v < m_n; ++v) {
    m_G.pts[v].r = randdblpos() * m_R;
    m_G.pts[v].phi = randdblpos() * 2.0 * M_PI;
  }
}

double SpringEmbedder::cooling() {
  // linear interpolation between 1.0 for i = 0 and 0.0 for i = steps
  return 1.0 - 1.0 / m_steps * m_i;

  // return sqrt((1.0 * m_steps - m_i) / m_steps);

  // return pow(m_i, -0.5);
}

double SpringEmbedder::forceStrength(double dist, bool adj) {
  // using hyperbolic probability with decreasing T for increasing i
  double T = 0.3 * m_R * max(0.05, cooling());  // 0.3
  // T = m_T;
  if (adj) {
    return (1 - Hyperbolic::prob(dist, m_R, T));
  } else {
    return (Hyperbolic::prob(dist, m_R, T));
  }
}

double SpringEmbedder::forceStrength1D(double delta_phi, double weight,
                                       double min_weight, double max_weight, int v1, int v2) {
  double opt_phi = (log(weight + 1) - log(min_weight + 1)) /
                   (log(max_weight + 1) - log(min_weight + 1)) * M_PI;
  opt_phi = min(M_PI, opt_phi);

  double diff = delta_phi - opt_phi;
  double sign = (0.0 < diff) - (diff < 0.0);
  diff = abs(diff);

  if (sign > 0 && diff > 0.5 * M_PI) {
    diff = M_PI - diff; // less influence for vertices on the opposite side
  }

  return sign * diff * diff;
}

void SpringEmbedder::computeForcesSlow(vector<HYPT>& forces) {
  // sum forces for all vertex pairs
  for (int v1 = 0; v1 < m_n; ++v1) {
    forces[v1].r = 0;
    forces[v1].phi = 0;
    for (int v2 = 0; v2 < m_n; ++v2) {
      if (v1 == v2) continue;

      HYPT f = force(m_G.pts[v1], m_G.pts[v2], m_adj_matrix[v1][v2]);

      forces[v1].r += f.r;
      forces[v1].phi += f.phi;
    }
  }
}

void SpringEmbedder::computeForces(vector<HYPT>& forces) {
  // use geometric data structure
  m_gds.init(&m_G);

  for (int v1 = 0; v1 < m_n; ++v1) {
    forces[v1].r = 0;
    forces[v1].phi = (0.8*cooling() + 0.2) * forces[v1].phi;
    // forces[v1].phi = 0;

    m_gds.forAllNodesDo(v1, [&](int v2, int count, HYPT pos, bool type) {
      (void)v2;
      (void)type;
      HYPT f = force(m_G.pts[v1], pos, false);
      forces[v1].r += count * f.r;
      forces[v1].phi += count * f.phi;
    });

    // compute attractive forces for neighbors of v1
    for (int v2 : m_G.edges[v1]) {
      // repair mistakes done before
      HYPT f = force(m_G.pts[v1], m_G.pts[v2], false);
      forces[v1].r -= f.r;
      forces[v1].phi -= f.phi;

      // compute the actual force
      f = force(m_G.pts[v1], m_G.pts[v2], true);
      forces[v1].r += f.r;
      forces[v1].phi += f.phi;
    }
  }
}

void SpringEmbedder::computeForces1D(vector<double>& forces,
                                     const vector<vector<double>>& weights,
                                     double min_weight, double max_weight) {
  for (int v1 = 0; v1 < m_n; ++v1) {
    forces[v1] = m_speed_factor *
                 ((1 - m_min_speed) * cooling() + m_min_speed) * forces[v1];
    // forces[v1] = (0.0*cooling() + 1) * forces[v1];
    for (int v2 = 0; v2 < m_n; ++v2) {
      if (v1 == v2) continue;

      forces[v1] += force1D(m_G.pts[v1].phi, m_G.pts[v2].phi, weights[v1][v2],
                            min_weight, max_weight, v1, v2);
    }
  }
}

void SpringEmbedder::computeAbsForces1D(vector<double>& forces,
                                        const vector<vector<double>>& weights,
                                        double min_weight, double max_weight) {
  for (int v1 = 0; v1 < m_n; ++v1) {
    forces[v1] = 0.0; // max(0.1, cooling()) * forces[v1];
    for (int v2 = 0; v2 < m_n; ++v2) {
      if (v1 == v2) continue;
      double f = force1D(m_G.pts[v1].phi, m_G.pts[v2].phi, weights[v1][v2],
                            min_weight, max_weight, v1, v2);
      forces[v1] += f < 0 ? -f : f;
    }
  }
}

void SpringEmbedder::minMaxWeights(const vector<vector<double>>& weights,
                                   double& min_weight, double& max_weight) {
  vector<double> tmp;
  max_weight = std::numeric_limits<double>::min();
  min_weight = std::numeric_limits<double>::max();
  for (auto vec : weights) {
    for (double weight : vec) {
      // get 2*median instead of the maximum
      tmp.push_back(weight);
      if (weight < min_weight) min_weight = weight;
    }
  }
  std::sort(tmp.begin(), tmp.end());
  max_weight = tmp[tmp.size()/2] * 2;
}


HYPT SpringEmbedder::force(const HYPT& p1, const HYPT& p2, bool adj) {
  HYPT f;

  double strength = forceStrength(HYPT::dist(p1, p2), adj);

  double tmp_phi = p2.phi - p1.phi;
  if (abs(tmp_phi) > M_PI) {
    tmp_phi += tmp_phi < 0 ? 2 * M_PI : -2 * M_PI;
  }

  double sign = (0.0 < tmp_phi) - (tmp_phi < 0.0);
  f.phi = strength * M_PI * sign * (adj ? 1 : -1);

  f.r = strength * m_R * (adj ? -1 : 1);

  return f;
}

double SpringEmbedder::force1D(double phi1, double phi2, double weight,
                               double min_weight, double max_weight, int v1, int v2) {
  double tmp_phi = phi2 - phi1;
  if (abs(tmp_phi) > M_PI) {
    tmp_phi += tmp_phi < 0 ? 2 * M_PI : -2 * M_PI;
  }

  double delta_phi = abs(tmp_phi);
  double strength = forceStrength1D(delta_phi, weight, min_weight, max_weight, v1, v2);

  double sign = (0.0 < tmp_phi) - (tmp_phi < 0.0);

  return strength * sign;
}

void SpringEmbedder::applyForces(const vector<HYPT>& forces, double max_step_r,
                                 double max_step_phi) {
  double max_r = 0;
  double max_phi = 0;
  for (int v = 0; v < m_n; ++v) {
    max_r = max(abs(forces[v].r), max_r);
    max_phi = max(abs(forces[v].phi), max_phi);
  }

  double factor_r = max_step_r / max_r;
  double factor_phi = max_step_phi / max_phi;

  for (int v = 0; v < m_n; ++v) {
    HYPT p = m_G.pts[v];

    p.r += factor_r * forces[v].r;
    p.phi += factor_phi * forces[v].phi;

    // r should be in [0, R]
    if (p.r < 0) {
      p.r = -p.r;
      p.phi += M_PI;
    }
    if (p.r > m_R) p.r = m_R;
    // phi should be in [0, 2pi]
    p.phi -= std::floor(p.phi / (2 * M_PI)) * 2 * M_PI;

    // set position
    m_gds.move(v, p);
  }
}

void SpringEmbedder::applyForces1D(const vector<double>& forces, double max_step) {
  double max_force = 0;
  for (int v = 0; v < m_n; ++v) {
    max_force = max(abs(forces[v]), max_force);
  }

  double factor = max_step/max_force;

  for (int v = 0; v < m_n; ++v) {
    double phi = m_G.pts[v].phi;
    
    phi += factor * forces[v];

    // phi should be in [0, 2pi]
    phi -= std::floor(phi / (2 * M_PI)) * 2 * M_PI;
    
    // set position
    // if (v == 3)
    m_G.pts[v].phi = phi;
  }
}




//  LocalWords:  img pdf
