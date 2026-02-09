#include "CustomDistributions.h"
#include <cmath>

// Normal
NormalDistribution::NormalDistribution(double xmin, double xmax, int npts,
                                       double mean, double stddev)
    : FiniteFunctions(xmin, xmax, npts), mu(mean), sigma(stddev) {}

double NormalDistribution::Evaluate(double x) const {
    return (1.0 / (sigma * sqrt(2 * M_PI))) *
           exp(-0.5 * pow((x - mu) / sigma, 2));
}

// Cauchy
CauchyDistribution::CauchyDistribution(double xmin, double xmax, int npts,
                                       double x0_, double gamma_)
    : FiniteFunctions(xmin, xmax, npts), x0(x0_), gamma(gamma_) {}

double CauchyDistribution::Evaluate(double x) const {
    return 1.0 / (M_PI * gamma *
           (1 + pow((x - x0) / gamma, 2)));
}
