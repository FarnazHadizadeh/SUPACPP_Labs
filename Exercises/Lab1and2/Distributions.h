#ifndef DISTRIBUTIONS_H
#define DISTRIBUTIONS_H

#include "FiniteFunctions.h"

/// Normal Gaussian distribution.
class NormalFunction : public FiniteFunctions {
private:
    double m_mu;
    double m_sigma;

public:
    NormalFunction(double xmin, double xmax,
                   double mu, double sigma,
                   std::size_t nPoints = 1000);

    double evaluate(double x) const override;
    std::string name() const override;
};

/// Cauchy–Lorentz distribution.
class CauchyLorentzFunction : public FiniteFunctions {
private:
    double m_x0;
    double m_gamma;

public:
    CauchyLorentzFunction(double xmin, double xmax,
                          double x0, double gamma,
                          std::size_t nPoints = 1000);

    double evaluate(double x) const override;
    std::string name() const override;
};

/// Crystal Ball distribution.
class CrystalBallFunction : public FiniteFunctions {
private:
    double m_xbar;
    double m_sigma;
    double m_alpha;
    double m_n;

    // Internal constants for normalisation.
    double m_N;
    double m_A;
    double m_B;

    void updateConstants();

public:
    CrystalBallFunction(double xmin, double xmax,
                        double xbar, double sigma,
                        double alpha, double n,
                        std::size_t nPoints = 1000);

    double evaluate(double x) const override;
    std::string name() const override;
};

#endif // DISTRIBUTIONS_H
