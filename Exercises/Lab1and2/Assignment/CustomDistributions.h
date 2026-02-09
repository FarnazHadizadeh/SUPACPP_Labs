#ifndef CUSTOMDISTRIBUTIONS_H
#define CUSTOMDISTRIBUTIONS_H

#include "FiniteFunctions.h"

// Normal Distribution
class NormalDistribution : public FiniteFunctions {
private:
    double mu, sigma;

public:
    NormalDistribution(double xmin, double xmax, int npts,
                       double mean, double stddev);
    double Evaluate(double x) const override;
};

// Cauchy Distribution
class CauchyDistribution : public FiniteFunctions {
private:
    double x0, gamma;

public:
    CauchyDistribution(double xmin, double xmax, int npts,
                       double x0_, double gamma_);
    double Evaluate(double x) const override;
};

#endif
