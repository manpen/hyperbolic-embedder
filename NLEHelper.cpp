#include "hyperboliclinear.hpp"
#include "NLEHelper.hpp"

#include <glog/logging.h>
#include <limits>

#include "plfit/powerlawCommon.h"
#include "random.hpp"

#define FOR(i,n) for(int (i)=0;(i)<(n);(i)++)
#define FORB(i,a,n) for(int (i)=(a);(i)<(n);(i)++)
#define FOREACH(it,c) for(__typeof((c).begin()) it=(c).begin();it!=(c).end();++it)
#define PB emplace_back
#define MP make_pair

void NLEHelper::estimateGlobalParameters() {

  vector<int> hist = H.degHisto();
  int missing_zeros = hist[1] - (hist[2] - hist[1]);

  n_orig = H.n + (missing_zeros > 0 ? missing_zeros : 0);
  m_orig = m;
  R = 2 * log(8 * n_orig * alpha * alpha * T
                  / (sin(M_PI * T) * (2. * m_orig / n_orig )
                     * (2*alpha-1) * (2*alpha-1)));
}

void NLEHelper::estimatePowerLaw() {
  plfit::VectorType degs;
  FOR (i,H.n)
    degs.PB(H.edges[i].size());
  plfit::VectorType results;
  plfit::Powerlaw::SingleFit(degs, results, false, false, 1.5, 0.01, 3);

  alpha = (results[0] - 1)/2;
  if (alpha <= 0.5) {
    LOG(WARNING) << "alpha estimated at " << alpha
                 << ", too low for embedding. Using alpha=0.505 instead.";
    alpha = 0.505;
  } else if (alpha >= 0.995) {
    LOG(WARNING) << "alpha estimated at " << alpha
                 << ", too high for embedding. Using alpha=0.995 instead.";
    alpha = 0.995;

  }
}

void NLEHelper::estimateHyperbolicParameters(const Graph& G, double *T,
                                             double *n, double *m,
                                             double *alpha, double *R,
                                             vector<double>* radial_coords) {
  // Setting T to a small arbitrary constant seems to produce good enough
  // results
  *T = 0.1;
  NLEHelper nle(G, *T);
  nle.radial_coords = radial_coords;

  nle.estimatePowerLaw();
  nle.estimateGlobalParameters();

  *n = nle.n_orig;
  *m = nle.m_orig;
  *alpha = nle.alpha;
  *R = nle.R;

  LOG(INFO) << "Estimated data:";
  LOG(INFO) << "alpha  = " << *alpha;
  LOG(INFO) << "T      = " << *T;
  LOG(INFO) << "n_orig = " << *n;
  LOG(INFO) << "m_orig = " << *m;
  LOG(INFO) << "R      = " << *R;

  if (radial_coords != nullptr)
    nle.estimateRadialCoordinates();
}

void NLEHelper::estimateRadialCoordinates() {
  FOR(i, H.n) {
    radial_coords->at(i) = 2 * log(2 * n_orig * alpha * T
                       / (sin(M_PI * T) * H.edges[i].size() * (alpha - 0.5)));
    if (radial_coords->at(i) > R)
      radial_coords->at(i) = R;
  }
}
