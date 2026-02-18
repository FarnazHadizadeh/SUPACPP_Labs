#include "FiniteFunctions.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <random>
#include <string>
#include <vector>

FiniteFunctions::FiniteFunctions(double xmin, double xmax, std::size_t nPoints)
    : m_xmin(xmin), m_xmax(xmax), m_nPoints(nPoints)
{
    if (m_xmin >= m_xmax) {
        std::swap(m_xmin, m_xmax);
    }
    if (m_nPoints == 0) {
        m_nPoints = 1000;
    }
}

double FiniteFunctions::xmin() const {
    return m_xmin;
}

double FiniteFunctions::xmax() const {
    return m_xmax;
}

std::size_t FiniteFunctions::nPoints() const {
    return m_nPoints;
}

double FiniteFunctions::evaluate(double x) const {
    // Default finite function: 1 / (1 + x^2)
    return 1.0 / (1.0 + x * x);
}

std::string FiniteFunctions::name() const {
    return "1 / (1 + x^2)";
}

double FiniteFunctions::integral(std::size_t nSamples) const {
    std::size_t n = (nSamples > 0 ? nSamples : m_nPoints);
    if (n == 0) {
        n = 1;
    }

    const double dx = (m_xmax - m_xmin) / static_cast<double>(n);
    double sum = 0.0;

    for (std::size_t i = 0; i <= n; ++i) {
        const double x = m_xmin + dx * static_cast<double>(i);
        const double w = (i == 0 || i == n) ? 0.5 : 1.0; // trapezoidal weights
        sum += w * evaluate(x);
    }

    return sum * dx;
}

double FiniteFunctions::normalisationFactor(std::size_t nSamples) const {
    const double I = integral(nSamples);
    if (I <= 0.0) {
        return 1.0;
    }
    return 1.0 / I;
}

void FiniteFunctions::writeFunctionToFile(const std::string &filename,
                                          bool normalised,
                                          std::size_t nSamples) const
{
    std::size_t n = (nSamples > 0 ? nSamples : m_nPoints);
    if (n == 0) {
        n = 1;
    }

    std::ofstream out(filename.c_str());
    if (!out) {
        return;
    }

    const double dx = (m_xmax - m_xmin) / static_cast<double>(n);
    const double norm = normalised ? normalisationFactor(nSamples) : 1.0;

    for (std::size_t i = 0; i <= n; ++i) {
        const double x = m_xmin + dx * static_cast<double>(i);
        const double y = evaluate(x) * norm;
        out << x << " " << y << "\n";
    }
}

std::vector<double> FiniteFunctions::metropolisSample(std::size_t nSamples,
                                                      double proposalSigma,
                                                      std::size_t burnIn,
                                                      std::mt19937 &rng,
                                                      bool normalise) const
{
    std::vector<double> samples;
    samples.reserve(nSamples);
    if (nSamples == 0) {
        return samples;
    }

    std::uniform_real_distribution<double> uniformX(m_xmin, m_xmax);
    std::uniform_real_distribution<double> uniform01(0.0, 1.0);
    std::normal_distribution<double> proposal(0.0, proposalSigma);

    double x = uniformX(rng); // starting point
    const double norm = normalise ? normalisationFactor() : 1.0;
    double fx = evaluate(x) * norm;
    if (fx < 0.0) {
        fx = 0.0;
    }

    const std::size_t totalSteps = burnIn + nSamples;

    for (std::size_t i = 0; i < totalSteps; ++i) {
        const double y = x + proposal(rng);

        if (y < m_xmin || y > m_xmax) {
            // Reject proposals that fall outside the defined range.
        } else {
            double fy = evaluate(y) * norm;
            if (fy < 0.0) {
                fy = 0.0;
            }

            double A = 1.0;
            if (fx > 0.0) {
                const double ratio = fy / fx;
                A = (ratio < 1.0) ? ratio : 1.0;
            }

            const double T = uniform01(rng);
            if (T < A) {
                x  = y;
                fx = fy;
            }
        }

        if (i >= burnIn) {
            samples.push_back(x);
        }
    }

    return samples;
}
