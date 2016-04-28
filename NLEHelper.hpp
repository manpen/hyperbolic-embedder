//
//  NLEHelper.hpp
//  Estimates global parameters and radii of nodes.
//

#pragma once

class NLEHelper {
 public:
  static void estimateHyperbolicParameters(const Graph& H, double *T,
                                           double* n, double* m, double* alpha,
                                           double* R,
                                           vector<double>* radial_coords);
 private:
  NLEHelper(const Graph& H, double T) : H(H), T(T) {
    m = 0;
    for (const auto& node : H.edges) m += node.size();
    m /= 2;
  }

  void estimatePowerLaw();
  void estimateGlobalParameters();
  void estimateRadialCoordinates();

  const Graph& H;
  double T;
  // Number of edges
  int m;

  // Estimated parameters
  double n_orig, m_orig;
  double alpha, R;

  vector<double>* radial_coords;
};
