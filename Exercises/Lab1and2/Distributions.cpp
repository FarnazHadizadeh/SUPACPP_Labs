#include "Distributions.h"

#include <cmath>
#include <string>

namespace {
    constexpr double pi = 3.14159265358979323846;
}

// ================= NormalFunction =================

NormalFunction::NormalFunction(double xmin, double xmax,
                               double mu, double sigma,
                               std::size_t nPoints)
    : FiniteFunctions(xmin, xmax, nPoints),
      m_mu(mu),
      m_sigma(sigma)
{
}

double NormalFunction::evaluate(double x) const {
    const double invSqrt2Pi = 1.0 / std::sqrt(2.0 * pi);
    const double z = (x - m_mu) / m_sigma;
    return invSqrt2Pi / m_sigma * std::exp(-0.5 * z * z);
}

std::string NormalFunction::name() const {
    return "Normal(mu=" + std::to_string(m_mu)
         + ", sigma=" + std::to_string(m_sigma) + ")";
}

// ================= CauchyLorentzFunction =================

CauchyLorentzFunction::CauchyLorentzFunction(double xmin, double xmax,
                                             double x0, double gamma,
                                             std::size_t nPoints)
    : FiniteFunctions(xmin, xmax, nPoints),
      m_x0(x0),
      m_gamma(gamma)
{
}

double CauchyLorentzFunction::evaluate(double x) const {
    const double t = (x - m_x0) / m_gamma;
    return 1.0 / (pi * m_gamma * (1.0 + t * t));
}

std::string CauchyLorentzFunction::name() const {
    return "Cauchy-Lorentz(x0=" + std::to_string(m_x0)
         + ", gamma=" + std::to_string(m_gamma) + ")";
}

// ================= CrystalBallFunction =================

CrystalBallFunction::CrystalBallFunction(double xmin, double xmax,
                                         double xbar, double sigma,
                                         double alpha, double n,
                                         std::size_t nPoints)
    : FiniteFunctions(xmin, xmax, nPoints),
      m_xbar(xbar),
      m_sigma(sigma),
      m_alpha(alpha),
      m_n(n),
      m_N(1.0),
      m_A(1.0),
      m_B(0.0)
{
    updateConstants();
}

void CrystalBallFunction::updateConstants() {
    const double absAlpha = std::fabs(m_alpha);
    const double expTerm  = std::exp(-0.5 * absAlpha * absAlpha);

    m_A = std::pow(m_n / absAlpha, m_n) * expTerm;
    m_B = m_n / absAlpha - absAlpha;

    const double C = (m_n / absAlpha) * (1.0 / (m_n - 1.0)) * expTerm;
    const double D = std::sqrt(pi / 2.0) *
                     (1.0 + std::erf(absAlpha / std::sqrt(2.0)));

    m_N = 1.0 / (m_sigma * (C + D));
}

double CrystalBallFunction::evaluate(double x) const {
    const double t = (x - m_xbar) / m_sigma;

    if (t > -m_alpha) {
        // Gaussian core
        return m_N * std::exp(-0.5 * t * t);
    } else {
        // Power-law tail
        return m_N * m_A * std::pow(m_B - t, -m_n);
    }
}

std::string CrystalBallFunction::name() const {
    return "CrystalBall(x0=" + std::to_string(m_xbar)
         + ", sigma=" + std::to_string(m_sigma)
         + ", alpha=" + std::to_string(m_alpha)
         + ", n=" + std::to_string(m_n) + ")";
}
