/**
 * @file
*
 * @author Manuel Penschuck
 * @copyright
 * Copyright (C) 2019 Manuel Penschuck
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * @copyright
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * @copyright
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
#include <algorithm>
#include <numeric>
#include <iostream>
#include <omp.h>
#include <cmath>
#include <cassert>
#include <algorithm>
#include <iostream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <sstream>

#include "graph.hpp"
#include "hyperbolic.hpp"
#include "hyperboliclinear.hpp"
#include "random.hpp"
#include "springembedder.hpp"
#include "NLEHelper.hpp"

#include <glog/logging.h>

#include "ScopedTimer.h"

#ifdef __unix__
#include "unistd.h"

static std::string hostname() {
    char tmp[128];
    if (gethostname(tmp, 128))
        return "n/a";
    return {tmp};
}

#else
static std::string hostname() {
    return "n/a";
}
#endif

/////////////////////////////////////////////////////////////
// copied from NetworKit
static double getExpectedDegree(double n, double alpha, double R) {
    double gamma = 2*alpha+1;
    double xi = (gamma-1)/(gamma-2);
    double firstSumTerm = exp(-R/2);
    double secondSumTerm = exp(-alpha*R)*(alpha*(R/2)*((M_PI/4)*pow((1/alpha),2)-(M_PI-1)*(1/alpha)+(M_PI-2))-1);
    double expectedDegree = (2/M_PI)*xi*xi*n*(firstSumTerm + secondSumTerm);
    return expectedDegree;
}

static double searchTargetRadiusForColdGraphs(double n, double k, double alpha, double epsilon) {
    double gamma = 2*alpha+1;
    double xiInv = ((gamma-2)/(gamma-1));
    double v = k * (M_PI/2)*xiInv*xiInv;
    double currentR = 2*log(n / v);
    double lowerBound = currentR/2;
    double upperBound = currentR*2;
    assert(getExpectedDegree(n, alpha, lowerBound) > k);
    assert(getExpectedDegree(n, alpha, upperBound) < k);
    do {
        currentR = (lowerBound + upperBound)/2;
        double currentK = getExpectedDegree(n, alpha, currentR);
        if (currentK < k) {
            upperBound = currentR;
        } else {
            lowerBound = currentR;
        }
    } while (abs(getExpectedDegree(n, alpha, currentR) - k) > epsilon );
    return currentR;
}

static double getTargetRadius(double n, double m, double alpha=1, double T=0, double epsilon = 0.01) {
    double result;
    double plexp = 2*alpha+1;
    double targetAvgDegree = (m/n)*2;
    double xiInv = ((plexp-2)/(plexp-1));
    if (T == 0) {
        double v = targetAvgDegree * (M_PI/2)*xiInv*xiInv;
        result = 2*log(n / v);
        result = searchTargetRadiusForColdGraphs(n, targetAvgDegree, alpha, epsilon);
    } else {
        double beta = 1/T;
        if (T < 1){//cold regime
            double Iinv = ((beta/M_PI)*sin(M_PI/beta));
            double v = (targetAvgDegree*Iinv)*(M_PI/2)*xiInv*xiInv;
            result = 2*log(n / v);
        } else {//hot regime
            double v = targetAvgDegree*(1-beta)*pow((M_PI/2), beta)*xiInv*xiInv;
            result = 2*log(n/v)/beta;
        }
    }
    return result;
}
/////////////////////////////////////////////////////////////

double benchmark(std::ostream& os, const std::string& host, unsigned iter, unsigned int n, unsigned int avgDeg, double alpha, double T, unsigned int seed = 0) {
    double time_total;
    long long num_edges;
    const double R = getTargetRadius(n, avgDeg / 2.0 * n, alpha, T);


    initRNG(seed);
    srand((unsigned int) seed);

    {
        ScopedTimer timer(time_total);
        auto H = HyperbolicLinear::linearSampling(n, R, alpha, T);
        num_edges = H.num_edges;
    }

    // Logging
    {
        std::stringstream ss;
        ss << "[CSV]"
           << host << ","
           << iter << ","
           << "Embedder,"
           << n << ","
           << avgDeg << ","
           << alpha << ","
           << T << ","
           << R << ","
           << time_total << ","
           << num_edges << ","
           << (2.0 * num_edges / n);

        os << ss.str() << std::endl;
    }

    return time_total;
}

int main(int argc, char* argv[]) {
    google::InitGoogleLogging(argv[0]);

    // Print Header
    std::cerr << "[CSV]"
       "host,"
       "iter,"
       "algo,"
       "n,"
       "avgDeg,"
       "alpha,"
       "T,"
       "R,"
       "TimeTotal,"
       "GenNumEdge,"
       "GenAvgDeg\n";

    unsigned seed = 0;
    const auto host = hostname();

    const unsigned n0 = 1e4;
    const unsigned nMax = 1e8;
    const unsigned steps_per_dec = 3;
    const double timeout = 100 * 1e3; // ms

    for(int iter = 0; iter < 5; iter++) {
        for (const double T : {0.0, 0.5, 0.9}) {
            for (double ple : {2.2, 3.0}) {
                const auto alpha = (ple - 1.0) / 2.0;
                unsigned int skip_n = nMax + 1;
                for (const auto avgDeg : {10, 100, 1000}) {

                    int ni = 0;
                    for (auto n = n0; n <= nMax; n = n0 * std::pow(10.0, 1.0 * ni / steps_per_dec), ++ni) {
                        std::cout << "\033[31miter=" << iter << ", PLE=" << ple << ", n=" << n << ", avgDeg=" << avgDeg << "\033[0m\n";

                        if (avgDeg * 20 > n) continue;

                        std::cout << "iter=" << iter << ", n=" << n << ", avgDeg=" << avgDeg << "\n";

                        double time;
                        if (n < skip_n) {
                            time = benchmark(std::cerr, host, iter, n, avgDeg, alpha, T, seed);
                            if (time > timeout) {
                                skip_n = n;
                                std::cout << " took too long\n";
                            }
                        } else {
                            std::cout << " skip_n = " << skip_n << "\n";
                        }
                        seed += 10;
                    }
                }
            }
        }
    }

    return 0;
}
