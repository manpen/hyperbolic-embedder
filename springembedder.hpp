#pragma once

#include <vector>
#include <string>
#include <iostream>

#include "geometricds.hpp"

using std::vector;
using std::string;
using std::function;
using std::ostream;
using std::istream;

class SpringEmbedder {
 public:
  SpringEmbedder(Hyperbolic& G);

  void run(int steps = 1000, int start = 1, int pictures = 1,
           bool rand_init_coords = true);

  // as the other run function, but considering only the first n
  // vertices
  void run(int n, int steps = 1000, int start = 1, int pictures = 1,
           bool rand_init_coords = true);

  // run the 1-dimensional spring embedder
  void run(const vector<vector<double>>& weights, int steps = 1000,
           int start = 1, int pictures = 1, bool rand_init_coords = true);

  // print all information necessary to later run the 1-dimensional
  // spring embedder
  void storeToStream(const vector<vector<double>>&weights, ostream& os);

  static void runFromStream(istream& is, vector<double>& res, double max_step,
                            double min_speed, double speed_factor);

  double score(const vector<vector<double>>& weights);

 private:

  // place vertices randomly
  void randomCoordinates();

  // strength of the force, depending on the distance between the points
  // and on whether the vertices are adjacent or not
  double forceStrength(double dist, bool adj);

  double forceStrength1D(double delta_phi, double weight, double min_weight,
                         double max_weight, int v1, int v2);

  // returns a factor between 1 and 0 depending on the current iteration
  // i and the total number of iterations
  double cooling();

  void computeForcesSlow(vector<HYPT>& forces);
  void computeForces(vector<HYPT>& forces);
  void computeForces1D(vector<double>& forces,
                       const vector<vector<double>>& weights, double min_weight,
                       double max_weight);
  void computeAbsForces1D(vector<double>& forces,
                          const vector<vector<double>>& weights,
                          double min_weight, double max_weight);

  HYPT force(const HYPT& p1, const HYPT& p2, bool adj);

  double force1D(double phi1, double phi2, double weight, double min_weight,
                 double max_weight, int v1, int v2);

  void applyForces(const vector<HYPT>& forces, double max_step_r,
                   double max_step_phi);

  void applyForces1D(const vector<double>& forces, double max_step);

  void minMaxWeights(const vector<vector<double>>& weights, double& min_weight, double& max_weight);

  Hyperbolic& m_G;
  int m_n;
  double m_T;

  int m_steps;
  int m_start;
  int m_i;
  int m_pictures;

  double m_max_step_1D = 0.55 * M_PI;
  double m_min_speed = 0.5;
  double m_speed_factor = 1.0;

  vector<vector<bool>> m_adj_matrix;

  double m_R;

  // GeometricDS<TreePoints> m_gds;
  GeometricDS<PointsFromCells> m_gds;
};

